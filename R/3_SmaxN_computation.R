#' Compute the Synchronised MaxN for a given second and a given camera
#'
#' This function computes the Synchronised MaxN (SmaxN) for a given second and
#' a given camera. It used the time dataframe gathering the time it takes for 
#' an individual of a given species to go from one camera to another. 
#' 
#' @param dist_df  a numerical dataframe containing the distance between each 
#'  pair of camera. There are as many rows as there are cameras and there are
#'  as many columns as there are cameras, thus the dataframe is symmetrical
#'  and the diagonal is filled with 0. \strong{Rows names and columns names
#'  must be cameras names}. \strong{BE CAREFUL that the cameras are 
#' the same and in the same order in the dist_df and in the abund_df!}
#' 
#' @param fish_speed a numerical value refering to the maximal speed of the 
#'  studied species. \strong{Speed must be given in meters per second}. If the
#'  computation of maxN values should not take into account fish speed (that is
#'  to say if the camera pooling is done at the second level), 
#'  \code{fish_speed = NULL}
#'  
#' @param abund_df a numerical dataframe containing the abundance of a 
#' given species across continuous time for several cameras. The columns refer
#' to the cameras and the rows refers to the time. \strong{Time must be given
#' in seconds and be continuous}. \strong{BE CAREFUL that the cameras are 
#' in the same order in the abund_df and the time_df!}.
#' 
#' @return if more than one camera: if the fish speed is taken into account, 
#' the function returns in a list: 
#' \strong{the Synchronised maxN (SmaxN) value} (the maximal value computed by 
#' gathering temporally synchronised cameras) , \strong{the maxN value} 
#' (the maximal value computed as the maximal abundance across all cameras and 
#' timesteps) , \strong{the maxN_cam values} (a vector of maximal values for 
#' each camera, thus the list contains as many elements as there are cameras) , 
#' \strong{the maxN_row values} (a vector of maximal values for each row by
#' taking the maximal value of cameras abundances on each row, thus the list 
#' contains as many elements as there are rows in the \code{abund_df} dataframe) 
#' and \strong{the SmaxN_row value} (the maximal value of the abundances summed
#' across cameras on each row). If fish speed is not taken into account then, 
#' the returned list does not contain the SmaxN value because the SmaxN 
#' value is the equal to the SmaxN_row value.
#' if only one camera: the function returns only \strong{the maxN value} 
#' (the maximal value computed as the maximal abundance across all timesteps)
#'  
#' @examples
#'  # Build distance dataframe for the example:
#'  dist_df_ex <- data.frame("A" = c(0, 2, 5, 5), "B" = c(2, 0, 5, 5), 
#'  "C" = c(5, 5, 0, 4), "D" = c(5, 5, 4, 0))
#'  rownames(dist_df_ex) <- c("A", "B", "C", "D")
#'  
#'  # Build distance dataframe for the example:
#'  abund_df_ex <- data.frame("A" = c(0, 1, 3, 7, 2, 2, 3, 0, 6, 2, 0, 1), 
#'                            "B" = c(2, 2, 2, 2, 0, 0, 0, 0, 1, 2, 1, 0), 
#'                            "C" = c(2, 0, 1, 0, 0, 4, 2, 2, 3, 0, 0, 4), 
#'                            "D" = c(0, 1, 0, 1, 0, 6, 1, 1, 6, 4, 2, 1))
#'  
#'  # Run the function:
#'  compute.max.abund(dist_df = dist_df_ex, 
#'                 fish_speed = 1.6, 
#'                 abund_df   = abund_df_ex)
#'                 
#' @export
#' 


compute.max.abund <- function(dist_df, fish_speed, abund_df) {
  
  
  # check that colnames of abund_df are the same as colnames and rownames ...
  # ... of dist_df:
  if (any(colnames(dist_df) != colnames(abund_df))) {
    stop("Cameras names are not the same  or not in the same order in the
         dist_df and the abund_df, please check.")
  }
  
  # check that colnames of abund_df are the same as colnames and rownames ...
  # ... of dist_df:
  if (any(colnames(dist_df) != rownames(dist_df))) {
    stop("Cameras names are not the same or not in the same order in the rows 
    names and the columns names of the dist_df, please check.")
  }
  
  # if there is more than  one camera:
  if (ncol(dist_df != 1)) {
    
    # if fish_speed is taken into account:
    if (! is.null(fish_speed)) {
      
      
      # get the time_df using the 1st function:
      time_df <- compute.cam.time(dist_df = dist_df, fish_speed = fish_speed)
      
      
      # get the length of the small interval build the lowest time between ...
      # ... camera pairs:
      small_UI <- min(apply(abund_df[, ], 1, function(x) min(x[x > 0])))
      # idem with max:
      big_UI <- max(abund_df)
      
      ## create two df for the SmaxN values for each timestep for the ...
      ## ... small_UI and big_UI:
      small_SmaxN_df <- as.data.frame(matrix(ncol = 2, nrow = 1))
      colnames(small_SmaxN_df) <- c("row", "SmaxN")
      big_SmaxN_df <- as.data.frame(matrix(ncol = 2, nrow = 1))
      colnames(big_SmaxN_df) <- c("row", "SmaxN")
      
      
      ## for each timestep, compute the SmaxN on small_UI and keep in memory ...
      ## ... the SmaxN for each timestep:
      ## idem for big_UI:
      
      for (i in (1:nrow(abund_df))) {
        
        # compute the SmaxN for the timestep and the small span:
        small <- compute.SmaxN.timestep(time_df = time_df,
                                        abund_df = abund_df,
                                        timestep = i,
                                        value = small_UI)

        # compute the SmaxN for the timestep and the big span:
        big <- compute.SmaxN.timestep(time_df = time_df,
                                        abund_df = abund_df,
                                        timestep = i,
                                        value = big_UI)
        
        # add values in the small and big df:
        small_SmaxN_df <- dplyr::add_row(row = rownames(abund_df)[i], SmaxN = small, small_SmaxN_df)
        big_SmaxN_df <- dplyr::add_row(row = rownames(abund_df)[i], SmaxN = big, big_SmaxN_df)
        
        # remove first rows filled with Na and use numerical format:
        small_SmaxN_df <- small_SmaxN_df[-1, ]
        big_SmaxN_df <- big_SmaxN_df[-1, ]
        
        small_SmaxN_df$SmaxN <- as.numeric(small_SmaxN_df$SmaxN)
        big_SmaxN_df$SmaxN <- as.numeric(big_SmaxN_df$SmaxN)
        
      } # end computation of Smaxn for each timestep and small and big spans
      
      
      # Compute the max of SmaxN of small spans for all timesteps:
      max_small <- max(small_SmaxN_df$SmaxN)
      
      
      
      ## Remove timesteps not to study ie the one with big_UI SmaxN < max_small:
      clean_big_SmaxN_df <- big_SmaxN_df[which(! big_SmaxN_df$SmaxN < max_small), ]
      # now this df contains the timesteps to study, in which the general SmaxN is
      
      # order the rows by decreasing order so that study intervals with the ...
      # ... biggest SmaxN first:
      order_big_SmaxN_df <- dplyr::arrange(clean_big_SmaxN_df, desc(SmaxN))
      
      
      # use the new function to compute SmaxN in a big interval 
      
      
      
      
      
      
      
    } # end if fish speed taken into account
    
  } # end if more than one camera
  
    
  # if there only one camera, return only maxN:
  if (ncol(dist_df) == 1) {
    return(list(maxN  = max(abund_df)))
  }

}
