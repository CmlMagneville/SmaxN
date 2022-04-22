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
#' @return the next value to be chosen that is to say the biggest value of the
#' bloc after the last value in max_bloc. If there are ex-aequos (ie several
#' cases of the bloc have the next highest value), the function returns a 
#' dataframe as max_bloc (Column names: value, cam_nm, timestep) with each row
#' being an ex-aequo value. If there is no ex-aequo, the function also returns a 
#' dataframe but with only one row.
#' 


search.next.value <- function(max_bloc, bloc) {
  
  # copy the bloc:
  bloc2 <- bloc
  
  # create an index which will take the values of the decreasing max values ...
  # ... of the bloc:
  value <- max(bloc2)
  
  # remove cameras which has already been chosen on the path:
  cam_used <- unique(max_bloc$cam_nm)
  bloc2 <- bloc2[, which(! colnames(bloc2) %in% cam_used)]

  # while the max value is already on the path, transform cases into NAs
  while(value %in% max_bloc$value) {
    bloc2[bloc2 == value] <- NA
    value <- max(bloc2, na.rm = TRUE)
  }
  
  
  # get the nb of times the highest value occurs:
  nb_occ <- sum(bloc == value)
  
  # create the df that will contain the coordinate(s) of the case(s) ...
  # ... which can be chosen because next highest value:
  next_value <- as.data.frame(matrix(ncol = 3, nrow = nb_occ))
  colnames(next_value) <- c("value", "cam_nm", "timestep")
    
  # index to fill the rows of this matrix 
  n <- 1
  
  # create a loop to get coordinates:
  for (i in (1:nrow(bloc))) {
    for (j in (1:ncol(bloc))) {

      if (bloc[i, j] == value) {
        next_value$value[n] <- value
        next_value$cam_nm[n] <- colnames(bloc)[j]
        next_value$timestep[n] <- rownames(bloc)[i]
        n <- n + 1
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
      ### PEUT POSER PROBLEME ICI EN FONCTION DE LA CLASSE DES TIME STEPS (CONVERTIR SI EN HMS)
      delta_real <- abs(as.numeric(max_bloc$timestep[i]) - as.numeric(next_possible2$timestep[j]) + 1)

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
#' @return 
#'
#' @example 
#' 
#'  abund_df <- data.frame("A" = c(9,8,3,3,3,3,3), 
#'  "B" = c(0,4,2,2,1,3,3), 
#'  "C" = c(0,0,0,0,1,1,1), 
#'  "D" = c(1,0,0,0,3,3,3))
#'  
#'  time_df <- data.frame("A" = c(0,5,6,6,7), "B" = c(5,0,3,3,3), 
#'  "C" = c(6,3,0,3,2), "D" = c(6,3,3,0,2), "E" = c(7, 3, 2, 2, 0))
#'  rownames(time_df) <- c("A", "B", "C", "D", "E")
#'  
#'  timestep <- 3
#'  
#'  value <- 6
#'  

compute.SmaxN.bigUI <- function(abund_df,
                                value,
                                timestep, 
                                time_df) {
  
  
  
  # get the dataframe which says in which order cameras should be looked at:
  cam_order_df <- compute.cam.order(time_df)
  
  
  # get the bloc (interval) to study:
  bloc <- abund_df[c(timestep:(value + timestep - 1)), ]
  
  # reduce to the number of cameras to use (if all 0):
  clean_bloc <- bloc[, colSums(bloc) != 0]
  
  # get the number of cameras:
  cam_nb <- ncol(clean_bloc)
  
  # if only one camera is kept:
  if (cam_nb == 1) {
    return(v = max(clean_bloc))
  }
  
  
  # if more than one camera to keep:
  if (cam_nb > 1) {
    
    
    # create a dataframe that will contain the max values chosen for each camera:
    max_bloc <- as.data.frame(matrix(ncol = 3, nrow = cam_nb))
    colnames(max_bloc) <- c("values", "cam_nm", "timestep")
    
    # create a counter to count camera number:
    n <- 1
    
    # while a max value not chosen for each camera:
    while (n <= cam_nb) {
      
     
      ## for the first camera: chose the highest value of the cleaned bloc 
      # ... NOTE: for now, if several highest values, chose the first one
      if (n == 1) {
        max <- max(clean_bloc)
        # chose the first one: EX AEQUO CHANGE
        max <- max[1]
        
        # create a loop to get coordinates:
        for (i in (1:nrow(clean_bloc))) {
          for (j in (1:ncol(clean_bloc))) {
            
            # span the cells to get the highest value: 
            if (clean_bloc[i, j] == max) {
              max_bloc$cam_nm[n] <- colnames(clean_bloc)[j]
              max_bloc$timestep[n] <- rownames(clean_bloc)[i]
              max_bloc$values[n] <- max
              break
            }
            
          }
          
          # once an occurrence has been found, stop:
          if (! is.na(max_bloc$timestep[n])) {
            break
          }
          
        }
        
        n <- n + 1
        
      } # end if n == 1 (chose the first value of the bloc)
      
      
      ## for the other cameras (chose values for cameras other than the first one):
      if (n > 1) {
        
        # search the new value (ex: 5):
        v <- search.next.value(max_bloc, clean_bloc)
        
        # know which cell is possible:
        possible <- next.possible(max_bloc[which(! is.na(max_bloc$values)), ], 
                                  next_possible = v, time_df)
      
          # if ex-aequo:
        if (length(v) > 1) { 
        
          # if at least one cell is possible, get the coord of the first one EX AEQUO:
          for (r in (1:nrow(possible))) {
            if (possible$possible[r] == TRUE) {
               max_bloc$values[n] <- possible$value[r]
               max_bloc$cam_nm[n] <- possible$cam_nm[r]
               max_bloc$timestep[n] <- possible$timestep[r]
               break
            }
          } # end get coord
          
         
        } # end if ex-aequo
        
        
        # if only one value:
        if (length(v) == 1) {
          max_bloc$values[n] <- possible$value[1]
          max_bloc$cam_nm[n] <- possible$cam_nm[1]
          max_bloc$timestep[n] <- possible$timestep[1]
        } 
        
        n <- n + 1
        
      } # end if  n > 1
      
      
    } # end while n < nb_cam
    
    
  } # end if more than one camera


  
  
  
}
  