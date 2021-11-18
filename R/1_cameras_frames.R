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
#'  and the diagonal is filled with 0.
#' 
#' @param fish_speed a numerical value refering to the mean speed of the 
#'  studied species. \strong{Speed must be given in meters per second}.
#' 
#' @return the function returns a dataframe containing the minimal time 
#'  needed for an individual of the studied species to go from a camera to 
#'  another camera.There are as many rows as there are cameras and there are
#'  as many columns as there are cameras, thus the dataframe is symmetrical
#'  and the diagonal is filled with 0.
#' 


compute.cam.time <- function(dist_df, fish_speed) {
  
}