################################################################################
##
## 3_Intermediate_functions.R
##
## This script gathers three functions used in the global SmaxN function 
## (script 4): 
## - one computing the order in which cameras should be classified,
## - one computing the frame of possible for a given timestep and a recursive
## - one computing all the possible paths given a timestep and the frame of 
##   possible and the biggest SmaxN found on these path
##
## 10/05/2022
##
## Camille Magneville
##
################################################################################




#' Compute the camera order
#'
#' This function rated the cameras in the following order: the first camera
#' is the central one (lowest mean distance across all other camera) (if several
#' camera possible, then classify in the order they come in the abund_df, it is 
#' not that important as we test all possible paths). Then other cameras are 
#' rated according to their minimal distances to the camers already chosen.
#' 
#' @param time_df a numerical dataframe containing the minimal time 
#'  needed for an individual of the studied species to go from a camera to 
#'  another camera.There are as many rows as there are cameras and there are
#'  as many columns as there are cameras, thus the dataframe is symmetrical
#'  and the diagonal is filled with 0. This dataframe is the output of the 
#'  \code{compute.cam.time} function.
#' 
#' @return a dataframe with one row containing the camera order to follow
#' 
#' @export
#' 
#' @examples 
#' time_df <- data.frame("A" = c(0,1,3,3), "B" = c(1,0,2,2), 
#'                       "C" = c(3,3,0,2), "D" = c(3,2,2,0))
#' rownames(time_df) <- c("A", "B", "C", "D")
#' 
#' order_cam <- cam.order(time_df)
#' 


cam.order <- function(time_df) {
  
  
  # create a df that contain the name of cameras already used in the right order:
  cam_df <- as.data.frame(matrix(ncol = ncol(time_df), nrow = 1))
  colnames(cam_df) <- c(1:ncol(time_df))
  
  # add a counter which count on which camera we are working
  # initialise with first camera:
  n <- 1
  
  # compute the number of cameras:
  nb_cam <- ncol(time_df)
  
  # create a time_df2 df which will be transformed thereafter:
  time_df2 <- time_df
  
  
  ## while all the cameras haven't been added:
  while(n <= nb_cam) {
    

    # if first camera to add, chose the central one:
    if (n == 1) {
      
      # change 0 -> NA so can compute mean with na.rm = TRUE:
      time_df2[time_df2 == 0] <- NA
      
      next_cam_df <- as.data.frame(apply(as.matrix(time_df2), 2, mean, na.rm = TRUE))
      next_cam_df <- tibble::rownames_to_column(next_cam_df, "cam")
      
      # chose the next camera to put in the next column as the one with minimal value:
      next_cam_nm <- next_cam_df$cam[which(next_cam_df[2] == min(next_cam_df[2]))]
      
      # get the number of cam to put:
      length_next_cam <- length(next_cam_nm)
      
      # fill the cam_df:
      cam_df[1, c(1:length_next_cam)] <- next_cam_nm
      
      # update n:
      n <- length(cam_df[, ! is.na(cam_df)]) + 1
      
      # update the time_df2 df with in row = camera left to add and column = added:
      time_df2 <- time_df2[! rownames(time_df2) %in% unlist(cam_df), colnames(time_df2) %in% unlist(cam_df)]
      time_df2 <- as.data.frame(time_df2)
      colnames(time_df2) <- cam_df[, which(! is.na(cam_df))]
      rownames(time_df2) <- rownames(time_df[which(! rownames(time_df) %in% unlist(cam_df)),])
      
    } # end if first camera to add
    
    
    # if not the first camera to add and not the last one:
    if (n > 1 & n < nb_cam) {
      
      # list of cameras to add:
      to_add <- setdiff(colnames(time_df), as.vector(unlist(cam_df[1, which(! is.na(cam_df))])))
      
      # get the minimal sum of the distance btw already used camera and camera to add:
      
      # if n = 2, then only one column in the time_df2 df so we don't have to compute the sum:
      if (n == 2) {
        min_value <- min(time_df2[to_add, ])
      }
      if (n > 2) {
        min_value <- min(apply(time_df2[to_add, ], 1, sum))
      }
      
      # loop on the rows of time_df2 to know which camera(s) has/have this min value:
      cam <- c()
      for (j in (1:nrow(time_df2))) {
        if (sum(time_df2[j, ]) == min_value) {
          cam <- append(cam, rownames(time_df2)[j])
        }
      } # end loop on rows to know which camera(s) to add next
      
      
      # fill the cam_df:
      cam_df[1, c(n:(n + length(cam) - 1))] <- cam
      
      # update n:
      n <- length(cam_df[, ! is.na(cam_df)]) + 1
      
      # update the time_df2 df with in row = camera left to add and column = added:
      time_df2 <- time_df[! rownames(time_df) %in% unlist(cam_df), colnames(time_df) %in% unlist(cam_df)]
      time_df2 <- as.data.frame(time_df2)

    } # end if not the first camera to add and not the last one
    
    
    # if the last camera to add:
    if (n == nb_cam) {
      
      # name of camera to add:
      to_add <- setdiff(colnames(time_df), as.vector(unlist(cam_df[1, which(! is.na(cam_df))])))
      
      # add the camera:
      cam_df[1, n] <- to_add
      
      # update the counter:
      n <- n + 1
      
    } # end if he last camera to add
    
    
  } ## end while all the camera haven't been added
  
  
  return("order_cam_df" <- cam_df)
  
  
}






#' Compute the frame of possible cells
#'
#' This function computes for a given timestep, the frame of possible cells:
#' the cells on which path can be computed. Practically the function adds NA
#' to the cells not on the frame of possibles in the abund_df. 
#' 
#' 
#' @param T a numerical value referring to the number of the row of the 
#' abundance dataframe studied (timestep)
#'
#' @param time_df a numerical dataframe containing the minimal time 
#'  needed for an individual of the studied species to go from a camera to 
#'  another camera.There are as many rows as there are cameras and there are
#'  as many columns as there are cameras, thus the dataframe is symmetrical
#'  and the diagonal is filled with 0. This dataframe is the output of the 
#'  \code{compute.cam.time} function.
#'  
#' @param abund_df a numerical dataframe containing the abundance of a 
#' given species across continuous time for several cameras. The columns refer
#' to the cameras and the rows refers to the time. \strong{Time must be given
#' in seconds and be continuous}. \strong{BE CAREFUL that the cameras are 
#' in the same order in the abund_df and the time_df!}. 
#' 
#' @return the abundance dataframe with NA on all cells except in the frame of 
#' possible around the studied timestep T.
#' 
#' @export
#' 
#' @examples 
#' 
#' abund_df <- data.frame("A" = c(9,8,3,3,3,3,3), "B" = c(0,4,2,2,1,3,3), 
#' "C" = c(0,0,0,0,1,1,1), "D" = c(1,0,0,0,3,3,3))
#' 
#' time_df <- data.frame("A" = c(0,4,4,3), "B" = c(4,0,2,2), 
#'                       "C" = c(4,2,0,2), "D" = c(3,2,2,0))
#' rownames(time_df) <- c("A", "B", "C", "D")
#' 
#' cam_order <- cam.order(time_df)
#' 
#' abund_df <- abund_df[, as.vector(unlist(cam_order))]
#' 
#' frame_possible <- frame.possible(T = 4, time_df, abund_df)
#' 


frame.possible <- function(T, time_df, abund_df) {
  
  
  # get the name of the first camera on which all paths will begin:
  cam1 <- colnames(abund_df)[1]
  
  # build a vector that contain the names of camera which are ok:
  # and initialise with the name of the first camera:
  cam_ok <- c(cam1)
  
  # initialise span object otherwise error when checking package:
  span <- 0
  span_start <- 0
  span_end <- 0
  
  # loop on all cameras, to complete each column of ...
  # ... the abund_df with NA where not possible frames:
  for (i in (1:ncol(abund_df))) {
    
    
    # for the first camera, only keep the studied timestep:
    if (i == 1) {
      
      rownames(abund_df) <- as.numeric(rownames(abund_df))
      abund_df[which(rownames(abund_df) > T), i] <- NA
      abund_df[which(rownames(abund_df) < T), i] <- NA
      
    }
    
    
    # if we study the second camera: what are its possible values?
    if (i == 2) {
      
      span <- time_df[cam1, colnames(abund_df)[i]]
      span_start <- T - span
      span_end <- T + span
      
      cam_ok <- append(cam_ok, colnames(abund_df)[i])
      
    } # end if second camera
    
    
    # if we don't study the second camera:
    if (i > 2) {
      
      # compute all the possible distances between cameras already updated ...
      # ... on the frame of possible and the studied camera:
      values <- c()
      
      # loop on camera already ok:
      for (j in (1:length(cam_ok))) {
        cam_done <- cam_ok[j]
        span <- time_df[cam_done, colnames(abund_df)[i]]
        values <- append(values, span)
      } # end loop to get distances btw studied cam i and cam already ok
      
      # get the span value which is:
      # min( min(dist_cam_i_all_cam_ok) , min(dist_cam_i_cam1))
      # = min(values):
      span <- min(values)
      span_start <- T - span
      span_end <- T + span
      
    } # end loop if don't study the second camera
    
    
    # check if span_start and span_end have problems (ie if at the beginning ...
    # ... or at the end of the abundance dataframe) and correct:
    
    # if pb at the beginning:
    if (T - span <= 0) {
      span_start <- 1
    }
   
    # if pb at the end:
    if (T + span > nrow(abund_df)) {
      span_end <- nrow(abund_df)
    }
    
    
    # Now let's "reduce" the abundance df by adding NAs outside the frame of ...
    # ... possible:
    
    if (span_start == 1 & span_end < nrow(abund_df)) {
      rownames(abund_df) <- as.numeric(rownames(abund_df))
      abund_df[which(rownames(abund_df) > span_end), i] <- NA
    }
    
    if (span_start > 1 & span_end == nrow(abund_df)) {
      rownames(abund_df) <- as.numeric(rownames(abund_df))
      abund_df[which(rownames(abund_df) < span_start), i] <- NA
    }
    
    if (span_start > 1 & span_end < nrow(abund_df)) {
      rownames(abund_df) <- as.numeric(rownames(abund_df))
      abund_df[which(rownames(abund_df) > span_end), i] <- NA
      abund_df[which(rownames(abund_df) < span_start), i] <- NA
    }

    cam_ok <- append(cam_ok, colnames(abund_df)[i])
    
  } # end loop on all cameras except the first one
  
  
  return("frame_possible" = abund_df)

}







#' Compute the biggest SmaxN for all the possible path of a given timestep
#'
#' This function computes 
#' 
#' 
#' @param T a numerical value referring to the number of the row of the 
#' abundance dataframe studied (timestep)
#'
#' @param time_df a numerical dataframe containing the minimal time 
#'  needed for an individual of the studied species to go from a camera to 
#'  another camera.There are as many rows as there are cameras and there are
#'  as many columns as there are cameras, thus the dataframe is symmetrical
#'  and the diagonal is filled with 0. This dataframe is the output of the 
#'  \code{compute.cam.time} function.
#'  
#' @param abund_df a numerical dataframe containing the abundance of a 
#' given species across continuous time for several cameras. The columns refer
#' to the cameras and the rows refers to the time. \strong{Time must be given
#' in seconds and be continuous}. \strong{BE CAREFUL that the cameras are 
#' in the same order in the abund_df and the time_df!}. 
#' 
#' @return the abundance dataframe with NA on all cells except in the frame of 
#' possible around the studied timestep T.
#' 
#' 
#' @examples 
#' 
#' abund_df <- data.frame("A" = c(9,8,3,3,3,3,3), "B" = c(0,4,2,2,1,3,3), 
#' "C" = c(0,0,0,0,1,1,1), "D" = c(1,0,0,0,3,3,3))
#' 
#' time_df <- data.frame("A" = c(0,4,4,3), "B" = c(4,0,2,2), 
#'                       "C" = c(4,2,0,2), "D" = c(3,2,2,0))
#' rownames(time_df) <- c("A", "B", "C", "D")
#' 
#' cam_order <- cam.order(time_df)
#' 
#' abund_df <- abund_df[, as.vector(unlist(cam_order))]
#' 
#' 
#' 
