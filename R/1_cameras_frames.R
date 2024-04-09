#' Compute frames on which SmaxN has to be computed
#'
#' This function computes a dataframe containing the minimal time 
#' needed for an individual of the studied species to go from a camera to 
#' another camera. It is based on the distance between each pair of cameras 
#' and the swim speed of the species.
#' 
#' @param dist_df a numerical dataframe containing the distance between each 
#'  pair of camera. There are as many rows as there are cameras and there are
#'  as many columns as there are cameras, thus the dataframe is symmetrical
#'  and the diagonal is filled with 0. \strong{Rows names and columns names
#'  must be cameras names}.
#' 
#' @param fish_speed a numerical value refering to the maximal speed of the 
#'  studied species. \strong{Speed must be given in meters per second}.
#' 
#' @return the function returns a dataframe containing the minimal time 
#'  needed for an individual of the studied species to go from a camera to 
#'  another camera.There are as many rows as there are cameras and there are
#'  as many columns as there are cameras, thus the dataframe is symmetrical
#'  and the diagonal is filled with 0.
#'  
#' @examples
#'  # Build distance dataframe for the example:
#'  dist_df_ex <- data.frame("A" = c(0, 2, 5, 5), "B" = c(2, 0, 5, 5), 
#'  "C" = c(5, 5, 0, 4), "D" = c(5, 5, 4, 0))
#'  rownames(dist_df_ex) <- c("A", "B", "C", "D")
#'  
#'  # Run the function:
#'  SmaxN::compute.cam.time(dist_df = dist_df_ex, fish_speed = 1.6)
#'  
#' @noRd


compute.cam.time <- function(dist_df, fish_speed) {
  
  
  # First convert cameras distance into a time given swim speed:
  time_df <- dist_df / fish_speed
  
  # And remove decimales:
  time_df <- as.data.frame(lapply(time_df, floor))
  
  # Give rownames:
  rownames(time_df) <- colnames(time_df)
  
  return(time_df)
}
