################################################################################
##
## 4_Pseudo_SmaxN_and_checks.R
##
## This script gathers ... functions used in the global SmaxN function 
## (script 5): 
## - 
## - 
## - 
##
## 12/05/2022
##
## Camille Magneville
##
################################################################################


#' Compute the max of cameras for a given second and a given span
#'
#' This function computes the maximal value of abundance for a given second and 
#' the "n" seconds after. "n" being set up by the "value" parameter .
#' 
#' @param time_df a numerical dataframe containing the minimal time 
#'  needed for an individual of the studied species to go from a camera to 
#'  another camera.There are as many rows as there are cameras and there are
#'  as many columns as there are cameras, thus the dataframe is symmetrical
#'  and the diagonal is filled with 0. This dataframe is the output of the 
#'  \code{compute.cam.time} function
#' 
#' @param abund_df a numerical dataframe containing the abundance of a 
#' given species across continuous time for several cameras. The columns refer
#' to the cameras and the rows refers to the time. \strong{Time must be given
#' in seconds and be continuous}. \strong{BE CAREFUL that the cameras are 
#' in the same order in the abund_df and the time_df!}.
#' 
#' @param timestep the number of the row to the studied timestep
#' 
#' @param value the value of the span of the interval on which the SmaxN is 
#' computed. Note: if the span is of 2 seconds, then the interval gathers 
#' 3 cases of the `abund_df`.
#' 
#' @return the function returns the SmaxN value for the given timestep and span
#'  
#' @examples
#'  # Build distance dataframe for the example:
#'  abund_df_ex <- data.frame("A" = c(0, 1, 3, 7, 2, 2, 3), 
#'                            "B" = c(2, 2, 2, 2, 0, 0, 0), 
#'                            "C" = c(2, 0, 1, 0, 0, 4, 2), 
#'                            "D" = c(0, 1, 0, 1, 0, 6, 1))
#'  
#'  # Run the function:
#'  
#' 


pseudoSmaxN.timestep <- function(time_df, abund_df, timestep, value) {
  
  
  
  # add rownames as a column so easier to manipulate:
  abund_df2 <- tibble::rownames_to_column(abund_df, "rownames")
  rownames(abund_df2) <- c(1:nrow(abund_df2))
  
  # if there is enough rows to build the "bloc" on which SmaxN should be ...
  # ... computed according to the "value" parameter:
  # ... note: number of cases = number of seconds + 1
  if (as.numeric(abund_df2[which(abund_df2$rownames == timestep), "rownames"]) + value - 1 <= nrow(abund_df2)) {
    
    # get the number of the row where the timestep is:
    row_nb <- as.numeric(abund_df2$rownames[which(abund_df2$rownames == timestep)])
    
    # create a vector that will contain one abundance value per camera to ...
    # ... sum them once the process is finished:
    max_cam_vect <- c()
    
    
    # span over the different cameras:
    for (j in (2:ncol(abund_df2))) {
      
      # get all the abundance values in the studied interval for the given cam:
      UI_abund_cam <- abund_df2[c(row_nb:(row_nb + value - 1)), j]
      
      # append the max value of this abundance data for the given cam and interval:
      max_cam_vect <- append(max(UI_abund_cam), max_cam_vect)
      
    } # end loop on cameras
    
    SmaxN <- sum(max_cam_vect)
    
  }
  
  
  # if there is not enough rows to build the "bloc" on which SmaxN should be ...
  # ... computed according to the "value" parameter:
  # ... note: number of cases = number of seconds + 1
  if ((nrow(abund_df2) - as.numeric(rownames(abund_df2[which(abund_df2$rownames == timestep), ])) + 1)
      < value + 1) {
    
    # get the number of the row where the timestep is:
    row_nb <- as.numeric(rownames(abund_df2[which(abund_df2$rownames == timestep), ]))
    
    # create a vector that will contain one abundance value per camera to ...
    # ... sum them once the process is finished:
    max_cam_vect <- c()
    
    
    # span over the different cameras:
    for (j in (2:ncol(abund_df2))) {
      
      # get all the abundance values in the studied interval for the given cam:
      UI_abund_cam <- abund_df2[c(row_nb:nrow(abund_df2)), j]
      
      # append the max value of this abundance data for the given cam and interval:
      max_cam_vect <- append(max(UI_abund_cam), max_cam_vect)
      
    } # end loop on cameras
    
    SmaxN <- sum(max_cam_vect)
    
  }
  
  return(SmaxN)
  
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
    for (i in (1:nrow(max_bloc[which(! is.na(max_bloc$value)), ]))) {
      
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
#' @examples 
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
  
    
    # loop on the cameras:
    for (j in (colnames(bloc))) {
      
      # if not the studied camera:
      if (! j %in% first_cell_df$cam_nm) {
        
        
        # get the coordinates of the studied cell:
        studied_cell_nm <- first_cell_df$cam_nm
        studied_cell_time <- as.numeric(first_cell_df$timestep)
        
        # get the span around the timestep of the cell for cam j according to ...
        # ... the distance to studied cell:
        # +1 because in time_df, the time_df shows seconds and here ...
        # ... we should have cells  (cell = seconds + 1)
        span <- time_df[first_cell_df$cam_nm, j] 
        
        # gather the values of cam j on which the max should be computed:
        
        start_span <- studied_cell_time - span
        end_span <- studied_cell_time + span
        
        # but we have to check that the timestep of the studied cell is not ...
        # ... on a border thus, timestep + span will not be possible:
        
        # if start pb:
        if (studied_cell_time - span <= rownames(bloc)[1]) { #equal because no 0 in the df rows
          
          # then the span must begin at the first row of the bloc:
          start_span <- rownames(bloc)[1]
          
        }
        
        # if end pb:
        if (studied_cell_time + span > rownames(bloc)[nrow(bloc)]) {
          
          # then the span must begin at the first row of the bloc:
          end_span <- rownames(bloc)[nrow(bloc)]
          
        }
        
        # max for the values of the cam j in the possible cells around the studied:
        max_cam <- max(bloc[which(rownames(bloc) %in% c(start_span:end_span)), j], na.rm = TRUE)
        vect_cam_values <- append(vect_cam_values, max_cam)
        
      } # end if not the studied cam
      
      
    } # end loop on cameras 
    
    # compute the sum of all the max value if start with the given cell:
    sum_cell <- sum(vect_cam_values) + as.numeric(first_cell_df$values)
    
    
    ## then compare this sum to the SmaxN_small_UI and decide if the cell should
    # ... be kept or not:
    if (sum_cell < SmaxN_small_UI) {
      first_cell_df$possible <- FALSE
    }  
    
    if (sum_cell >= SmaxN_small_UI) {
      first_cell_df$possible <- TRUE
    }
    
    
  
  # return the updated first_cell_df with TRUE or FALSE for each cell:
  return(first_cell_df)
  
}


