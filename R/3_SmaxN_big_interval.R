################################################################################
##
## 3_SmaxN_big_interval.R
##
## This script gathers three functions used to find the path to the biggest 
## SmaxN in a given bloc (interval of time begining at a given timestep).
##
## 19/04/2022
##
## Camille Magneville
##
################################################################################



#' Give the next value on a given bloc given the first highest values.
#'
#' This function computes next number value to be chosen in a time bloc. It just gives 
#' the next highest value, another function informs if this value can be chosen
#' or not given the time it takes a fish to go between camera pairs.
#' 
#' @param max_bloc a dataframe containing in rows the number already chosen on
#' the path and their coordinates ie. on which column (camera name) and row 
#' (time) they are. Values are ranged fom the biggest (first rows) to the 
#' lowest (last rows). Column names: value, cam_nm, timestep
#'  
#' @param bloc a dataframe being the abundance dataframe restricted to the
#' studied bloc (ie a timestep + span). 
#' 
#' @param n the number of the column on which the max must be computed
#' 
#' @return the next value to be chosen that is to say the biggest value of the
#' bloc after the last value in max_bloc. If there are ex-aequos (ie several
#' cases of the bloc have the next highest value), the function returns a 
#' dataframe as max_bloc (Column names: value, cam_nm, timestep) with each row
#' being an ex-aequo value. If there is no ex-aequo, the function also returns a 
#' dataframe but with only one row.
#' 


search.next.value <- function(max_bloc, bloc, n) {
  
  # copy the bloc:
  bloc2 <- bloc
  
  # create an index which will take the values of the decreasing max values ...
  # ... of the bloc:
  # take the max of the column:
  value <- max(bloc2[, n], na.rm = TRUE)
  
  # # remove cameras which has already been chosen on the path:
  # cam_used <- unique(max_bloc$cam_nm[which(! is.na(max_bloc$values))])
  # bloc2 <- bloc2[, which(! colnames(bloc2) %in% cam_used)]
  
  # get the nb of times the highest value occurs:
  nb_occ <- sum(bloc[, n] == value, na.rm = TRUE)
  
  # create the df that will contain the coordinate(s) of the case(s) ...
  # ... which can be chosen because next highest value:
  next_value <- as.data.frame(matrix(ncol = 3, nrow = nb_occ))
  colnames(next_value) <- c("values", "cam_nm", "timestep")
  
  # index to fill the next_value df:
  k <- 1
  
  # create a loop to get coordinates:
  for (i in (1:nrow(bloc))) {
    
    if (! is.na(bloc[i, n])) {
      
      if (bloc[i, n] == value) {
        next_value$values[k] <- value
        next_value$cam_nm[k] <- colnames(bloc)[n]
        next_value$timestep[k] <- i               # check that should be i or rownames(bloc)[i]
        k <- k + 1
      }
    }
  }
  
  # return:
  return(next_value)
  
}



#' Tell if the next value is possible or not given time between cameras
#'
#' This function a TRUE/FALSE value if the next value is possible. If there
#' are ex-aequo values, begin with the first one in the \code{next_possible} 
#' dataframe, use the other until it works.
#' 
#' @param max_bloc a dataframe containing in rows the number already chosen on
#' the path and their coordinates ie. on which column (camera name) and row 
#' (time) they are. Values are ranged fom the biggest (first rows) to the 
#' lowest (last rows). Column names: value, cam_nm, timestep
#' 
#' @param next_possible  a dataframe with column names: value, cam_nm, timestep 
#' with a row for the information about the next possible value. If several cells 
#' have the same value (ex-aequo), the dataframe has one row for each cell. This
#' dataframe is obtained through the \code{search.next.value} function.
#' 
#' @param time_df a numerical dataframe containing the minimal time 
#'  needed for an individual of the studied species to go from a camera to 
#'  another camera.There are as many rows as there are cameras and there are
#'  as many columns as there are cameras, thus the dataframe is symmetrical
#'  and the diagonal is filled with 0. This dataframe is the output of the 
#'  \code{compute.cam.time} function.
#' 
#' @return the next value to be chosen that is to say the biggest value of the
#' bloc after the last value in max_bloc which is possible according to the 
#' distance between cameras. Note: For one possible value, check if possible
#' with already on the path values and stop when one is not possible. This
#' function tests for all possible values (all ex-aequo) so we have a df 
#' showing TRUE or FALSE for all ex-aequo given the values already on the path.
#' 

next.possible <- function(max_bloc, next_possible, time_df) {
  
  
  # add a column to the next_possible df which contain either TRUE or FALSE ...
  # ... if the case can be chosent according to the distance between camera:
  next_possible2 <- next_possible
  next_possible2$possible <- rep(NA, nrow(next_possible2))
  
  # loop on the possible values:
  for (j in (1:nrow(next_possible2))) {
    
    # create a vector gathering if possible or not
    possible_vect <- c()
    
    # loop on the already chosen max: those already on the "path":
    for (i in (1:nrow(max_bloc))) {
      
      # name of the camera of the value on the path:
      cam_path <- max_bloc$cam_nm[i]
      cam_next <- next_possible2$cam_nm[j]
      
      # if the cameras are different:
      if (cam_path != cam_next) {
        
        # compute the time delta which is ok (+1 because one timestep is two cases):
        delta_ok <- time_df[cam_path, cam_next] + 1
      }
      
      
      # if the cameras the same: normalement pas possible car fct qui cherche valeur verifie que cam differentes
      if (cam_path == cam_next) {
        
        delta_ok <- 0
        
      }
      
      # compute the real timespan between the value on the path and the next possible:
      # add +1 because time_df is in seconds so number of timesteps + 1
      ### PEUT POSER PROBLEME ICI EN FONCTION DE LA CLASSE DES TIME STEPS (CONVERTIR SI EN HMS)
      delta_real <- abs(as.numeric(max_bloc$timestep[i]) - as.numeric(next_possible2$timestep[j])) + 1

      # check the deltas:
      if (delta_real > delta_ok) {
        possible <- FALSE
        
        # add the value to the vector gathering if possible or not
        possible_vect <- append(possible_vect, possible)
        
        # stop the loop: no need to check for other cameras because at least not possible for one:
        break
      }
      
      if (delta_real <= delta_ok) {
        possible <- TRUE
        # add the value to the vector gathering if possible or not
        possible_vect <- append(possible_vect, possible)
      }
      
    } # end loop on the already chosen max
    
    # if all cameras are ok, chose this value and stop the scanning of all possible values ...
    # ... here chose the first which work but to change in case of ax-aequo:
    if (! FALSE %in% possible_vect) {
      next_possible2$possible[j] <- TRUE
    }
    else {
      next_possible2$possible[j] <- FALSE
    }
    
    
  } # end loop on possible values
  
  return(next_possible2)
  
}



#' Tell if possible cells for the first cameras are possibly leading to the 
#' SmaxN or if they are not.
#'
#' This function completes the \code{first_cell_df} which contains the cells
#' for a given value and the first camera by adding TRUE or FALSE if the 
#' given cell can lead to the SmaxN or will not lead to the SmaxN by comparing
#' to the maximal SmaxN computed on small UI.
#' 
#' @param time_df a numerical dataframe containing the minimal time 
#'  needed for an individual of the studied species to go from a camera to 
#'  another camera.There are as many rows as there are cameras and there are
#'  as many columns as there are cameras, thus the dataframe is symmetrical
#'  and the diagonal is filled with 0. This dataframe is the output of the 
#'  \code{compute.cam.time} function.
#'  
#' @param first_cell_df a dataframe with value, cam_nm and timesteps of 
#' possible cells for the first camera
#'
#' @param SmaxN_small_UI a numeric value referring to the maximal value of the 
#' SmaxN for the small UI
#'
#' @param bloc the studied big UI 
#' 
#' @return the completed \code{first_cell_df}
#'
#' @example 
#' 
#'  abund_df <- data.frame("A" = c(9,8,3,3,3,3,3), 
#'  "B" = c(0,4,2,2,1,3,3), 
#'  "C" = c(0,0,0,0,1,1,1), 
#'  "D" = c(1,0,0,0,3,3,3))
#'  
#'  time_df <- data.frame("A" = c(0,5,6,6), "B" = c(5,0,3,3), 
#'  "C" = c(6,3, 0, 3), "D" = c(6,3,3,0))
#'  rownames(time_df) <- c("A", "B", "C", "D")
#'  


first.cam.possible <- function(time_df, first_cell_df,
                               SmaxN_small_UI, bloc) {
  
  
  # create the vector:
  vect_cam_values <- c()
  
  # add a possible value to the first_cell_df:
  first_cell_df$possible <- rep(NA, nrow(first_cell_df))

  # loop on all cells:
  for (i in (1:nrow(first_cell_df))) {
    
    
    ## first get for each cameras the possible values and store them in a vector:
    
    # loop on the cameras:
    for (j in (colnames(bloc))) {
      
      # if not the studied camera:
      if (! j %in% first_cell_df$cam_nm) {
        
        
        # get the coordinates of the studied cell:
        studied_cell_nm <- first_cell_df$cam_nm[i]
        studied_cell_time <- as.numeric(first_cell_df$timestep[i])
        
        # get the span around the timestep of the cell for cam j according to ...
        # ... the distance to studied cell:
        # +1 because in time_df, the time_df shows seconds and here ...
        # ... we should have cells  (cell = seconds + 1)
        span <- time_df[first_cell_df$cam_nm[i], j] 
        
        # gather the values of cam j on which the max should be computed:
        
        start_span <- studied_cell_time - span
        end_span <- studied_cell_time + span

        # but we have to check that the timestep of the studied cell is not ...
        # ... on a border thus, timestep + span will not be possible:
        
        # if start pb:
        if (studied_cell_time - span <= 0) { #equal because no 0 in the df rows
          
          # then the span must begin at the first row of the bloc:
          start_span <- 1
          
        }
        
        # if end pb:
        if (studied_cell_time + span > nrow(bloc)) {
          
          # then the span must begin at the first row of the bloc:
          end_span <- nrow(bloc)
          
        }
        
        # max for the values of the cam j in the possible cells around the studied:
        max_cam <- max(bloc[start_span:end_span, j])
        vect_cam_values <- append(vect_cam_values, max_cam)
        
      } # end if not the studied cam
      
      
    } # end loop on cameras 
    
    # compute the sum of all the max value if start with the given cell:
    sum_cell <- sum(vect_cam_values) + as.numeric(unique(first_cell_df$values))
      
      
    ## then compare this sum to the SmaxN_small_UI and decide if the cell should
    # ... be kept or not:
    if (sum_cell < SmaxN_small_UI) {
      first_cell_df$possible[i] <- FALSE
    }  
    
    if (sum_cell >= SmaxN_small_UI) {
      first_cell_df$possible[i] <- TRUE
    }
    

  } # end loop on all cells 
  
  
  # return the updated first_cell_df with TRUE or FALSE for each cell:
  return(first_cell_df)
  
}



#' Compute the order in which cameras must be looked at in a bloc (smallest
#' to highest distances)
#'
#' This function computes a dataframe which contains one or several rows which...
#' camera names in the order they must be looked at. If several camera ...
#' pairs have the same distance, then their will be several rows... 
#' If all camera pairs have different distances, then there will be one row ... 
#' with only one camera order to test
#' 
#' @param time_df a numerical dataframe containing the minimal time 
#'  needed for an individual of the studied species to go from a camera to 
#'  another camera.There are as many rows as there are cameras and there are
#'  as many columns as there are cameras, thus the dataframe is symmetrical
#'  and the diagonal is filled with 0. This dataframe is the output of the 
#'  \code{compute.cam.time} function.
#' 
#' @return 
#'
#' @example 
#' 
#'  abund_df <- data.frame("A" = c(9,8,3,3,3,3,3), 
#'  "B" = c(0,4,2,2,1,3,3), 
#'  "C" = c(0,0,0,0,1,1,1), 
#'  "D" = c(1,0,0,0,3,3,3))
#'  
#'  time_df <- data.frame("A" = c(0,5,6,6), "B" = c(5,0,3,3), 
#'  "C" = c(6,3, 0, 3), "D" = c(6,3,3,0))
#'  rownames(time_df) <- c("A", "B", "C", "D")
#'  



compute.cam.order <- function(time_df) {
  
  # create a df that contain the name of cameras already used in the order list:
  cam_df <- as.data.frame(matrix(ncol = ncol(time_df), nrow = 1))
  colnames(cam_df) <- c(1:ncol(time_df))
  
  time_df2 <- time_df
  
  # counter to count the number of cameras:
  n <- 1
  
  # while each camera has not been studied:
  
  while(n <= ncol(time_df)) {
    
    # create a df with only the cameras which have not been selected already:
    time_df2 <- time_df[! rownames(time_df) %in% unlist(cam_df), ! colnames(time_df) %in% unlist(cam_df)]
    
    # if the last camera to put in the df
    if (is.numeric(time_df2)) {
      
      # get the camera name:
      cam_done <- unique(unlist(cam_df))
      time_df3 <- time_df
      time_df3 <- tibble::rownames_to_column(time_df3, "cam")
      cam_toadd <- time_df3$cam[! colnames(time_df) %in% cam_done]
      
      # add the camera:
      cam_df[, n] <- rep(cam_toadd, nrow(cam_df))
      
      break
    }
    
    # change 0 -> NA so can compute mean with na.rm = TRUE:
    time_df2[time_df2 == 0] <- NA
    
    next_cam_df <- as.data.frame(apply(as.matrix(time_df2), 2, mean, na.rm = TRUE))
    next_cam_df <- tibble::rownames_to_column(next_cam_df, "cam")
    
    # chose the next camera to put in the next column as the one with minimal value:
    next_cam_nm <- next_cam_df$cam[which(next_cam_df[2] == min(next_cam_df[2]))]
    
    
    # if only one camera with the lowest time & if there is one list to test so far:
    if (length(next_cam_nm) == 1 & (nrow(cam_df) == 1)) {
      
      cam_df[1, n] <- next_cam_nm
      n <- n + 1
      
    }
    
    # if only one camera with the lowest time & if there is more than one list to test so far:
    if (length(next_cam_nm) == 1 & (nrow(cam_df) > 1)) {
      
      for (j in (1:(length(next_cam_nm)-1))) {
        
        cam_df <- dplyr::add_row(cam_df, cam_df[1, ])
        cam_df[, n] <- next_cam_nm
        n <- n + 1
        
        
      }
    }
    
    
    # if two cameras with the lowest time:
    if (length(next_cam_nm) > 1) {
      
      
      # compute all combination of the camera which have the same distance:
      comb <- combinat::permn(next_cam_nm)
      
      # create as many rows as there are similar distances:
      j <- 1
      
      # add enough rows to be able to add all the possible combinations:
      cam_df <- cam_df[rep(c(1:nrow(cam_df)), length(comb)),]
      # reset row names:
      rownames(cam_df) <- c(1:nrow(cam_df))
      
      # add combinations to the df:
      for (m in (1:nrow(cam_df))) {
        cam_df[m, c(n:(n + length(next_cam_nm) - 1))] <- comb[[j]]
        j <- j + 1
      }
      
      
      n <- n + length(next_cam_nm)
    }
  }
  
  return(cam_df)
  
}

 


#' Compute the SmaxN for a timespan ie here the big interval of a given timestep
#'
#' This function computes the SmaxN value of the big interval. It choses the 
#' highest value for each camera taking into account the distance between the
#' cameras.
#' 
#' @param abund_df a numerical dataframe containing the abundance of a 
#' given species across continuous time for several cameras. The columns refer
#' to the cameras and the rows refers to the time. \strong{Time must be given
#' in seconds and be continuous}. \strong{BE CAREFUL that the cameras are 
#' in the same order in the abund_df and the time_df!}. 
#' 
#' @param value the span of the interval on which the SmaxN should
#' be computed (given in cell number) ie if in the `time_df` the number is 2
#' then value is 3 (cells).
#' 
#' @param timestep the timestep on which the interval begins
#' 
#' @param time_df a numerical dataframe containing the minimal time 
#'  needed for an individual of the studied species to go from a camera to 
#'  another camera.There are as many rows as there are cameras and there are
#'  as many columns as there are cameras, thus the dataframe is symmetrical
#'  and the diagonal is filled with 0. This dataframe is the output of the 
#'  \code{compute.cam.time} function.
#'  
#' @param SmaxN_small_UI a numeric value rrefering to the span of the interval
#'   to study.
#' 
#' @return 
#'
#' @example 
#'  " define the abundance df:
#'  abund_df <- data.frame("A" = c(9,8,3,3,3,3,3), 
#'  "B" = c(0,4,2,2,1,3,3), 
#'  "C" = c(0,0,0,0,1,1,1), 
#'  "D" = c(1,0,0,0,3,3,3))
#'  # define the time_df
#'  time_df <- data.frame("A" = c(0,5,6,6), "B" = c(5,0,3,3), 
#'  "C" = c(6,3,0,3), "D" = c(6,3,3,0))
#'  rownames(time_df) <- c("A", "B", "C", "D")
#'  # define the timestep:
#'  timestep <- 3
#'  # define the interval's span:
#'  value <- 4
#'  

compute.SmaxN.bigUI <- function(abund_df,
                                value,
                                timestep, 
                                time_df,
                                SmaxN_small_UI) {
  
  # create the vector with the SmaxN values for each terminated path:
  SmaxN_vect <- c()
  
  # get the dataframe which says in which order cameras should be looked at:
  cam_order_df <- compute.cam.order(time_df)
  
  # create a vector which contains SmaxN values computed for different paths...
  # ... and/or different bloc: what is of interest is its value not ...
  # ... its position:
  SmaxN_value <- c()
  
  # get the bloc (interval) to study:
  bloc <- abund_df[c(timestep:(value + timestep - 1)), ]
  
  
  # LOOP on all the camera order that should be studied:
  for (m in (1:nrow(cam_order_df))) {
    
    print(paste0("camera order", sep = " ", m, sep = " ", "on", sep = " ", nrow(cam_order_df)))
    
    # get the camera order:
    cam_order <- cam_order_df[m, ]
    
    # order the bloc based on camera names:
    bloc <- bloc[, as.vector(unlist(cam_order))]
    
    # set rownames:
    rownames(bloc) <- c(1:nrow(bloc))
    
    # reduce to the number of cameras to use (if all 0):
    clean_bloc <- bloc[, colSums(bloc) != 0]
    
    # get the number of cameras:
    cam_nb <- ncol(clean_bloc)
    
    # if no camera to keep:
    if (cam_nb == 0) {
      return("max_SmaxN_timestep" = 0)
    }
    
    # if only one camera is kept:
    if (cam_nb == 1) {
      return("max_SmaxN_timestep" = max(clean_bloc))
    }
    
    
    # if more than one camera to keep:
    if (cam_nb > 1) {
      
      
      # create a dataframe that will contain the max values chosen for each camera:
      path <- as.data.frame(matrix(ncol = 3, nrow = 1))
      colnames(path) <- c("values", "cam_nm", "timestep")
      # timestep here is the row number in the studied bloc
      
      # create a counter to count camera number done on the path:
      n <- 1
      
      print(paste0("find max for camera", sep = " ", n))
      
      
      # create a vector that contain values already tested to begin the path
      # (if a value can not work because of distances between cameras):
      first_val_vect <- c()
      
      
      # create vectors that contain values already tested for a given camera (n):
      for (j in (2:cam_nb)) {
        # create a vector that contain the already tested value for the given camera:
        assign(paste0("already_tested_max_n", sep = "_", j), c())
      }
      
      # LOOP while a max value not chosen for each camera:
      while (n <= cam_nb) {
        
        
        ## for the first camera: 
        if (n == 1) {
          
          
          # create a list for a tried path containing for each camera (column) ...
          # ... the df of possible cells (all TRUE for distance):
          list <- list()
          
          
          # get all the values of this camera to be able to stop the process if all
          # ... possible cells of this camera in the bloc have been checked:
          abund_first_cam <- unique(bloc[, n])
          
          # if all values of the first cam have been checked and do not work or 
          # ... do not allow to build a path with a good SmaxN, return ...
          # ... NA:
          if (length(first_val_vect) == length(abund_first_cam)) {
            return("max_SmaxN_timestep" = NA)
          }
          
          
          # create a df with cells with highest values
          # either one row if one cell or several rows if several cells:
          # build it with as many rows than there are cells in the bloc (to ...
          # ... many rows and will delete thereafter):
          first_cell_df <- as.data.frame(matrix(ncol = 3, 
                                                nrow = nrow(clean_bloc)*ncol(clean_bloc)))
          colnames(first_cell_df) <- c("values", "cam_nm", "timestep")
          
          
          # LOOP while the max for the bloc and the first camera (first column ...
          # ... is already in the values which have been tested, find another one:
          max <- max(clean_bloc[, 1], na.rm = TRUE)
          
          
          while(max %in% first_val_vect) {
            
            # search a new max of the column to test:
            # create a new df and replace values which have been ...
            # ... already tested by NA (so can be several cells even if ...
            # ... one value because several cells can have the same value in ...
            # ... the studied column):
            interm_bloc <- clean_bloc
            interm_bloc[interm_bloc[, 1] %in% first_val_vect, 1] <- NA
            
            # search a new max:
            max <- max(interm_bloc[, 1], na.rm = TRUE)
            
          } # end of while max is in already tested values:
          
          
          # LOOP to get the coordinates of cell(s) which have this max value:
          for (i in (1:nrow(clean_bloc))) {
            
            # span the cells to get the highest value: 
            if (clean_bloc[i, 1] == max) {
              
              cam_nm <- colnames(clean_bloc)[1]
              timestep <- rownames(clean_bloc)[i]
              values <- max
              first_cell_df <- dplyr::add_row(first_cell_df, cam_nm = cam_nm,
                                              timestep = timestep,
                                              values = values)
              
            }
          } # end loop to get coord of cell(s) having the max value
          
          # remove NA rows:
          first_cell_df <- first_cell_df[which(! is.na(first_cell_df$values)), ]
          
          print(paste0("cells for camera", sep = " ", n, sep = " ", "and max =", sep = " ", max,
                       sep = " ", "is:")) 
          print(first_cell_df)
          
          # check which first cell can be kept and which cell can not ...
          # ... lead to the SmaxN and thus remove these cells thereafter:
          first_cell_df2 <- first.cam.possible(time_df, first_cell_df, SmaxN_small_UI, bloc)
          
          print(paste0("max for camera", sep = " ", n, sep = " ", "is", sep = " ", unique(first_cell_df2$values)))
          
          # if at least one cell is possible:
          if (any(TRUE %in% first_cell_df2$possible)) {
            
            # only keep possible values and remove possible value:
            first_cell_df <- first_cell_df2[which(first_cell_df2$possible == TRUE), ]
            first_cell_df <- first_cell_df[, -ncol(first_cell_df)]
            
            # take the first row of first_cell_df in the path df:
            path <- dplyr::add_row(path, first_cell_df[1, ])
            
            # remove the first row which is filled with NAs:
            path <- path[-1, ]
            
            print("path so far is:")
            print(path)
            
            # update the list gathering all cells for values on the path:
            list[[n]] <- first_cell_df
            names(list[n]) <- n
            
            # update the counter to search the next camera: 
            n <- n + 1
            
          }
          
          
          # if no cell is possible, add the value to values tested and not possible:
          if (! any(TRUE %in% first_cell_df2$possible)) {
            
            # update the values already tested for the (n) cam
            first_val_vect <- append(first_val_vect, unique(first_cell_df2$values))
            
            print(paste0("cam", sep = " ", unique(first_cell_df2$cam_nm), sep = " ",
                         "and value", sep = " ", unique(first_cell_df2$values), sep = " ",
                         "do not work. Abandon path."))
            
            # n does not evolve and is still = 1
            
          }
          
          
        } # end if n == 1 (chose the first value of the bloc)
        
        
        
        
        ## for the other cameras (chose values for cameras other than the first one):
        if (n > 1) {
          
          print(paste0("find max for camera", sep = " ", n))
          
          # search the new value (ex: 5):
          v <- search.next.value(max_bloc = path[which(! is.na(path$values)), ],
                                 bloc = clean_bloc,
                                 n = n)
          
          # # while the max has already been tested:
          while(unique(v$value) %in% get(paste0("already_tested_max_n", sep = "_", n))) {
            #
            #   # search a new max of the column to test:
            #   # create a new df and replace values which have been ...
            #   # ... already tested by NA (so can be several cells even if ...
            #   # ... one value because several cells can have the same value in ...
            #   # ... the studied column):
            interm_bloc <- clean_bloc
            interm_bloc[which(interm_bloc[, n] == unique(v$value)), n] <- NA
            
            #   # search a new max:
            v <- search.next.value(max_bloc = path[which(! is.na(path$values)), ],
                                   bloc = interm_bloc,
                                   n = n)
            
          } # end of while max is in already tested values:
          
          
          # know which cell is possible:
          possible <- next.possible(path[which(! is.na(path$values)), ], 
                                    next_possible = v, time_df)
          
          
          print(paste0("possible max cells for camera", sep = " ", n, sep = " ", "and max =", sep = " ", max,
                       sep = " ", "are:"))
          print(possible)
          
          
          # if some cell(s) is/are possible:
          if (any(TRUE %in% possible$possible)) {
            
            # remove cells that are not possible:
            possible <- possible[which(possible$possible[] == TRUE), ]
            
            
            # compute the sum: already chosen values + max of columns which ...
            # ... have not yet been chosen to compare with the highest ...
            # ... SmaxN of the small UI:
            
            # if still several cameras left:
            if (n < cam_nb) {
              
              if (is.data.frame(clean_bloc[, (n + 1):ncol(clean_bloc)])) {
                
                S <- sum(path$values, na.rm = TRUE) + unique(v$value) +
                  sum(apply(clean_bloc[, (n + 1):ncol(clean_bloc)], 2, max, na.rm= TRUE),
                      na.rm = TRUE)
              }
              
            }
            
            # if only one cameras left to add:
            if (n < cam_nb) {
              
              if (is.numeric(clean_bloc[, (n + 1):ncol(clean_bloc)])) {
                
                S <- sum(path$values, na.rm = TRUE) + unique(v$value) +
                  max(clean_bloc[, (n + 1):ncol(clean_bloc)], na.rm= TRUE)
                
              }
              
            }
            
            # if study the last camera:
            if (n == cam_nb) {
              
              S <- sum(path$values, na.rm = TRUE) + unique(v$value) 
              
            }
            
            
            # test if S is ok for this path with the new value or if ...
            # ... S < small:
            if (S >= SmaxN_small_UI) {
              sum <- TRUE
            }
            
            if (S < SmaxN_small_UI) {
              sum <- FALSE
            }
            
            print(paste0("sum value is", sep = " ", sum))
            
            # if the sum == TRUE: then we can take this new value 
            # and we take the first cell having this value. (If does not work ...
            # ... after we will take the second cell if there is one etc):
            if (sum == TRUE) {
              
              path$timestep <- as.numeric(path$timestep)
              path <- dplyr::add_row(path, possible[1, -ncol(possible)])
              
              print("path so far is:")
              print(path)
              
              # search the next max value:
              n <- n + 1
              
              
              # update the list gathering all cells for values on the path:
              list[[n - 1]] <- possible[, -ncol(possible)]
              names(list[n - 1]) <- n - 1
              
            }
            
            
            # if the sum == FALSE: then taking this value will not lead ...
            # ... to a "good" SmaxN: it will be lowest than SmaxN_small_UI: ...
            # ... we have a problem on the step before (the value we have chosen ...
            # ... do not lead to a path with a sufficiently high SmaxN):
            # ... so we must abandon the path (as we start by using max of ...
            # ... each column) and start a new one with n <- 1:
            
            if (sum == FALSE) {
              
              n <- 1
              
              # create a clean path df:
              path <- as.data.frame(matrix(ncol = 3, nrow = 1))
              colnames(path) <- c("values", "cam_nm", "timestep")
              
            } # end if sum == FALSE
            
            
          } # end of "if some cells are possible"
          
          
          # if no cell is possible, chose another value for previous camera:
          if (! any(TRUE %in% possible$possible)) {
            
            print("if no cell possible for cam n")
            
            
            # if it is the second value (camera), if we remove the first row than no dataframe so put NAs:
            if ((n-1) == 1)  {
              
              print("if second cam")
              
              
              # if the value of the 1st camera has several cell, we take the next
              # cell:
              if (nrow(list[[n-1]]) > 1) {
                
                # remove the first cell of the list because does not work:
                list[[n-1]] <- list[[n-1]][- 1, ]
                # take the second cell of the list which is now the first one: 
                path[n - 1, ] <- list[[n-1]][1, ]
                
                print("path so far is:")
                print(path)
                
              } # end if the previous value has several cells
              
              
              # if the previous cam and value have only one cell (or if we have ...
              # ... already removed all possible with previous loop), ...
              # ... then we will not use the n value:
              if (nrow(list[[n-1]]) == 1) {
                
                print("if the 1 cam has only one cell")
                
                # recreate a unfilled path df:
                path <- as.data.frame(matrix(ncol = 3, nrow = 1))
                colnames(path) <- c("values", "cam_nm", "timestep")
                
                # add the value to the first_cam_vect because value for the ...
                # ... first camera does not work: tested all 1st cam cells ...
                # ... and all 2nd cam cells:
                first_val_vect <- append(first_val_vect, unique(list[[n-1]]$values))
                
                # check a new value for the first cam:
                n <- n - 1
                
                print(paste0("n is equal to", sep = " ", n))
                
              }

            }
            
            # if it is not the second camera:
            if ((n-1) > 1) {
              
              print("if not second cam")
              print(length(list))
              print(list[[n-1]]) # pb still several values
              
              # if the previous value has only one cell (or if we have ...
              # ... already removed all possible with previous loop), ...
              # ... then we will not use the n value:
              if (nrow(list[[n-1]]) == 1) {

                print("if the n-1 cam has only one cell")

                # update the values already tested for the (n) cam
                assign(paste0("already_tested_max_n", sep = "_", n), append(
                  get(paste0("already_tested_max_n", sep = "_", n)),
                  unique(v$values)))

                print(paste0("n is equal to", sep = " ", n))

              }
              
              # if the value of the (n-1) camera has several value, we take the next
              # cell:
              if (nrow(list[[n-1]]) > 1) {
                
                # remove the first cell of the list because does not work:
                list[[n-1]] <- list[[n-1]][- 1, ]
                # take the second cell of the list which is now the first one: 
                path[n - 1, ] <- list[[n-1]][1, ]
                
                print("path so far is:")
                print(path)

              } # end if the previous value has several cells
              
            } #  end if not the second camera and no cell left on the first
            
          } # end if no cell is possible
          
          
        } # end if  n > 1
        
        
      } # end while n < nb_cam
      
      
    } # end if more than one camera
    
    
    # update the vector with the SmaxN values for each terminated path:
    SmaxN_vect <- append(SmaxN_vect, sum(path$values, na.rm = TRUE))
    print(SmaxN_vect)
    
  } # end loop on camera order (build one bloc for each order to follow)
  
  return("max_SmaxN_timestep" = max(SmaxN_vect))

}