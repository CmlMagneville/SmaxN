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
#' with a row for the information about the next possible value. If several cases 
#' have the same value (ex-aequo), the dataframe has one row for each case. This
#' dataframe is obtained through the \code{search.next.value function.
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
      delta_real <- abs(max_bloc$timestep[i] - as.numeric(next_possible2$timestep[j]) + 1)

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
