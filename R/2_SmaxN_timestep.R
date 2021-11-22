#' Compute the Synchronised MaxN for a given second and a given camera
#'
#' This function computes the Synchronised MaxN (SmaxN) for a given second and
#' a given camera. It used the time dataframe gathering the time it takes for 
#' an individual of a given species to go from one camera to another. 
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
#' @param cam_nm the name of the camera for which the SmaxN must be computed
#' 
#' @param timestep the number of the row to the studied timestep
#' 
#' @return the function returns a list of the Synchronised MaxN (SmaxN) values
#' for all the cameras and a given second (timestep).
#'  
#'  @examples
#'  # Build distance dataframe for the example:
#'  abund_df_ex <- data.frame("A" = c(0, 1, 3, 7, 2, 2, 3), 
#'                            "B" = c(2, 2, 2, 2, 0, 0, 0), 
#'                            "C" = c(2, 0, 1, 0, 0, 4, 2), 
#'                            "D" = c(0, 1, 0, 1, 0, 6, 1))
#'  
#'  # Run the function:
#'  
#' 


compute.SmaxN.timestep <- function(time_df, abund_df, cam_nm, timestep) {
  
  
  vect_max <- c()
  
  # check that abund_df and time_df have columns ordered likely:
  
  
  # Do a loop to browse all the cameras (ie all the columns of abund_df):
  for(i in (1:ncol(abund_df))) {
    
    
    # if the i th camera is the one studied colnames(abund_df)[i] == cam_nm ...
    # ... then the frame on this column is equal to the minimal time for the ...
    # ... cam_nm row:
    if (colnames(abund_df)[i] == cam_nm) {
      
      start_frame <- timestep - min(time_df[cam_nm, 
                                            which(time_df[cam_nm, ] != 0)])
      stop_frame <- timestep + min(time_df[cam_nm, 
                                           which(time_df[cam_nm, ] != 0)])
      max <- max(abund_df[c(start_frame:stop_frame), i])
      vect_max <- append(vect_max, max)
    }
    
    # else:
    if (colnames(abund_df)[i] != cam_nm) {
      
      # extract max of each column in the progression frame:
      start_frame <- timestep - time_df[cam_nm, colnames(abund_df)[i]]
      stop_frame <- timestep + time_df[cam_nm, colnames(abund_df)[i]]
      max <- max(abund_df[c(start_frame:stop_frame), i])
      vect_max <- append(vect_max, max)
    }
    
  }
  
  
  # compute the sum of the max retrieved for each camera:
  sum_max_cam <- 0
  for (j in vect_max) {
    sum_max_cam <- sum_max_cam + j
  }
    
  return(sum_max_cam)
  
}
