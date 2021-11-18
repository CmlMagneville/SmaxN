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
#' for all the cameras and a given second.
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
  
  
  # Do a loop to browse all the cameras (ie all the columns of abund_df):
  for(i in (1:ncol(abund_df))) {
    
    # extract max of each column in the progression frame:
    max(abund_df[c(timestep), i])
    
  }
  
  
  
}
