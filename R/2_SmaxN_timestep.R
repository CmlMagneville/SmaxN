#' Compute the Synchronised MaxN for a given second and a given span
#'
#' This function computes the Synchronised MaxN (SmaxN) for a given second and 
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
#' @return the function returns a list of the Synchronised MaxN (SmaxN) values
#' for all the cameras and a given second (timestep).
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


compute.SmaxN.timestep <- function(time_df, abund_df, timestep, value) {
  
  
  
  # add rownames as a column so easier to manipulate:
  abund_df2 <- tibble::rownames_to_column(abund_df, "rownames")
  rownames(abund_df2) <- c(1:nrow(abund_df2))
  
  # if there is enough rows to build the "bloc" on which SmaxN should be ...
  # ... computed according to the "value" parameter:
  # ... note: number of cases = number of seconds + 1
  if ((nrow(abund_df2) - as.numeric(rownames(abund_df2[which(abund_df2$rownames == timestep), ])) + 1)
       >= value + 1) {
    
    # get the number of the row where the timestep is:
    row_nb <- as.numeric(rownames(abund_df2[which(abund_df2$rownames == timestep), ]))
    
    # create a vector that will contain one abundance value per camera to ...
    # ... sum them once the process is finished:
    max_cam_vect <- c()
    
    
    # span over the different cameras:
    for (j in (2:ncol(abund_df2))) {
      
      # get all the abundance values in the studied interval for the given cam:
      UI_abund_cam <- abund_df2[c(row_nb:(row_nb + value)), j]
      
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
