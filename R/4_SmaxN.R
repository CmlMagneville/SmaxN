#' Compute the Synchronised MaxN (SmaxN) (highest value)
#'
#' This function computes the Synchronised MaxN (SmaxN) for a set on n cameras
#' with know distances between them and for a given species which maximal speed
#' is known.
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
#' @param paral.option a logical value referring to the 
#' parallelisation process: if \code{FALSE}, then no parallelisation is 
#' realised, if \code{TRUE}, then parallelisation is done on 
#' (number of cores - 1) (so ok to do anything else on your machine). 
#' Currently, parallelisation is only possible on Windows computers.
#' 
#' @return if more than one camera: if the fish speed is taken into account, 
#' the function returns in a list: 
#' \strong{the Synchronised maxN (SmaxN) value} (the maximal value computed by 
#' gathering temporally synchronised cameras) , \strong{the maxN value} 
#' (the maximal value computed as the maximal abundance across all cameras and 
#' timesteps) , \strong{the maxN_cam values} (a vector of maximal values for 
#' each camera, thus the list contains as many elements as there are cameras) , 
#' \strong{the maxN_tiemstep values} (a vector of maximal values for each row by
#' taking the maximal value of cameras abundances on each row, thus the list 
#' contains as many elements as there are rows in the \code{abund_df} dataframe) 
#' and \strong{the SmaxN_timestep value} (the maximal value of the abundances summed
#' across cameras on each row). If fish speed is not taken into account then, 
#' the returned list does not contain the SmaxN value because the SmaxN 
#' value is the equal to the SmaxN_row value.
#' if only one camera: the function returns only \strong{the maxN value} 
#' (the maximal value computed as the maximal abundance across all timesteps)
#'  
#'  
#'  
#' @examples
#'  # Build distance dataframe for the example:
#'  dist_df_ex <- data.frame("A" = c(0, 2, 5, 5), "B" = c(2, 0, 5, 5), 
#'  "C" = c(5, 5, 0, 4), "D" = c(5, 5, 4, 0))
#'  rownames(dist_df_ex) <- c("A", "B", "C", "D")
#'  # Build distance dataframe for the example:
#'  abund_df_ex <- data.frame("A" = c(0, 1, 3, 7, 2, 2, 3, 0, 6, 2, 0, 1), 
#'                            "B" = c(2, 2, 2, 2, 0, 0, 0, 0, 8, 2, 1, 0), 
#'                            "C" = c(2, 0, 1, 0, 0, 4, 2, 2, 3, 0, 0, 4), 
#'                            "D" = c(0, 1, 0, 1, 0, 6, 1, 1, 6, 4, 2, 1))
#'  
#'                 
#' @export
#' 



SmaxN.computation <- function(abund_df, fish_speed, dist_df, paral.option) {
  
  
  #### Checks and basic manips
  
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
  
  # put rownames of the abundance_df as number:
  rownames(abund_df) <- c(1:nrow(abund_df)) 
  
  
  #### Parallelisation def
  
  if (paral.option == TRUE) {
    
    # returns the number of available hardware threads:
    nb_cores <- parallel::detectCores(logical = TRUE) 
    
    #  allocate this number of available cores to R and register:
    cl <- parallel::makeCluster(nb_cores - 1)  
    doParallel::registerDoParallel(cl) 
    
  }
  
 
  #### Begin the SmaxN serach
  
  ### Checks if no easy solution (SmaxN_small_UI):
  
  # if there is more than  one camera:
  if (ncol(dist_df) != 1) {
    
    # if fish_speed is taken into account:
    if (! is.null(fish_speed)) {
      
      
      # get the time_df using the 1st function:
      time_df <- compute.cam.time(dist_df = dist_df, fish_speed = fish_speed)
      
      
      # get the length of the small interval build the lowest time between ...
      # ... camera pairs:
      small_UI <- min(apply(time_df[, ], 1, function(x) min(x[x > 0])))
      
      ## create one df for the false "SmaxN" values for each timestep for the ...
      ## ... small_UI and big_UI:
      small_SmaxN_df <- as.data.frame(matrix(ncol = 2, nrow = 1))
      colnames(small_SmaxN_df) <- c("row", "SmaxN")
      

      ## for each timestep, compute the SmaxN on small_UI and keep in memory ...
      ## ... the SmaxN for each timestep:

      for (i in (1:nrow(abund_df))) {
        
        # compute the SmaxN for the timestep and the small span:
        small <- pseudoSmaxN.timestep(time_df = time_df,
                                        abund_df = abund_df,
                                        timestep = i,
                                        value = small_UI)
        
        # add values in the small and big df:
        small_SmaxN_df <- dplyr::add_row(row = rownames(abund_df)[i], SmaxN = small, small_SmaxN_df)
        
        small_SmaxN_df$SmaxN <- as.numeric(small_SmaxN_df$SmaxN)

      } # end computation of Smaxn for each timestep and small and big spans
      
      # remove first rows filled with Na and use numerical format:
      small_SmaxN_df <- small_SmaxN_df[-1, ]

      # Compute the max of SmaxN of small spans for all timesteps:
      max_small <- max(small_SmaxN_df$SmaxN)
      
      # order the rows by decreasing order so that study intervals with the ...
      # ... biggest SmaxN first:
      order_SmaxN_df <- dplyr::arrange(small_SmaxN_df, desc(SmaxN))
      
      
    #### order the abundance df so cameras ok:
      
      # compute the order in whch cameras should be viewed:
      cam_order <- cam.order(time_df)
      
      # order the abundance columns:
      abund_df2 <- abund_df
      abund_df2 <- abund_df2[, as.vector(unlist(cam_order))]
      abund_df <- abund_df2
      
    ### Search the highest SmaxN for each timestep to study:
      
      
      # create a vector that will contain the SmaxN of each big UI:
      SmaxN_vect <- c(max_small)
      
      # create a vector that contain one path to have the SmaxN:
      path_saved <- c()
      
      b <- 1
      
      while (b <= nrow(order_SmaxN_df)) {
        
        print(paste0("!!!!!! Starting row ",
                     sep = " ", b, sep = " ", "on", 
                     sep = " ", nrow(order_SmaxN_df)))
        print(paste0((b/nrow(order_SmaxN_df))*100, sep = "", "%"))
        
        
        
        # compute the frame of possible around the studied timestep:
        frame_possible_df <- frame.possible(T = as.numeric(order_SmaxN_df$row[b]), 
                                            time_df, abund_df)
        
        # reduce the frame_possible_df by deleting rows with all NAs:
        frame_possible_df <- frame_possible_df[rowSums(is.na(frame_possible_df)) != ncol(frame_possible_df), ]
        
        # remove columns with only 0 (but not first camera if the cell = 0):
        #diff_0 <- frame_possible_df[, colSums(frame_possible_df, na.rm = TRUE) != 0]
        #frame_possible_df <- cbind(frame_possible_df[, 1], diff_0)
        
        

        # check that with the first cell of this timestep we can have a "good"
        # SmaxN:
        
        first_cam_nm <- colnames(frame_possible_df)[1]
        
        first_cell_df <- data.frame("values" = abund_df[as.numeric(order_SmaxN_df$row[b]), 
                                                        first_cam_nm],
                                    "cam_nm" = first_cam_nm,
                                    timestep = as.numeric(order_SmaxN_df$row[b]))
        
        
        # check if first cell can be kept or not ...
        # ... lead to the SmaxN and thus remove these cells thereafter:
        first_cell_df2 <- first.cam.possible(time_df, first_cell_df, 
                                             SmaxN_small_UI = max_small, 
                                             bloc = frame_possible_df)
        
        # if the first cell is not possible (does not lead to "good" SmaxN):
        # then stop the process for this timestep
        if (first_cell_df2$possible == FALSE) {
          b <- b + 1
        }
        
        
        #### Create the path df:
        
        if (first_cell_df2$possible == TRUE) {
          
          # create the path_df:
          path_df <-  as.data.frame(matrix(ncol = 3, nrow = ncol(abund_df)))
          colnames(path_df) <- c("value", "cam_nm", "timestep")
          
          # add the first cam:
          path_df[1, ] <- c(first_cell_df2$values, 
                            first_cell_df2$cam_nm,
                            first_cell_df2$timestep)
          
          
          #### Search the first path to work on 
          # NOTE: not the one with the highest values
          # possible on each camera so we can get rid of a lot of other rows
          # thereafter: because algorithm takes too long to find! As first vers...:
          
          # Thus chose the row of the timestep:
          for (j in (2:nrow(path_df))) {
            
            cam_nm <- colnames(frame_possible_df)[j]
            timestep <- first_cell_df2$timestep
            value <- frame_possible_df[timestep, cam_nm]
            
            # complete the path_df
            path_df[j, ] <- c(value, cam_nm, timestep)
            
          }
          
          # Compute the sum for this first path and compare to SmaxN_small_UI:
          path_df$value <- as.numeric(path_df$value)
          S <- sum(path_df$value)
          
          if (S > max_small) {
            max_small <- S
          }
          
          # remove the last camera of the path,thus initialising the
          # recursive fct and running backward:
          path_df[nrow(path_df), ] <- rep(NA, 3)
          
          #### Compute all possible paths and keep in mind the highest SmaxN:
          # try without parallelisation
          
          # compute the SmaxN of the big interval of the given timestep:
          T <- rownames(frame_possible_df[!is.na(frame_possible_df[, 1]), ])
          T <- as.numeric(T)
          
          print("Starting recursive ")
          
          value <- recursive.paths(T = T, 
                               frame_possible_df, 
                               n = ncol(frame_possible_df), 
                               path_df = path_df, 
                               SmaxN_small_UI = max_small,
                               time_df)
          v <- value$"SmaxN_timestep"
          
          print("Ending recursive ")
          
          # add the SmaxN of the given timestep to the SmaxN vect:
          SmaxN_vect <- append(SmaxN_vect, v)
          path_saved <- append(path_saved, unlist(value$"one_path"))
          
          b <- b + 1
          
        } # end if first cell possibly leads to a "good" SmaxN

      } # end while
      
      
      # Compute the general SmaxN:
      SmaxN <- max(SmaxN_vect)
      
      # Compute other metrics:
      maxN <- max(abund_df)
      maxN_cam <- apply(abund_df, 2, max)
      maxN_timestep <- apply(abund_df, 1, max)
      
      vect_maxN_sum <- apply(abund_df, 1, sum)
      SmaxN_timestep <- max(vect_maxN_sum) 
      
      
      # return:
      return_list <- list("maxN" = maxN, 
                          "SmaxN" = SmaxN,
                          "path_saved" = path_saved,
                          "SmaxN_timestep" = SmaxN_timestep,
                          "maxN_cam" = maxN_cam,
                          "maxN_timestep" = maxN_timestep)
      
      return(return_list)
  
    } # if there is a fish speed
    
    
    if (is.null(fish_speed)) {
      
      # Compute other metrics:
      maxN <- max(abund_df)
      maxN_cam <- apply(abund_df, 2, max)
      maxN_timestep <- apply(abund_df, 1, max)
      
      vect_maxN_sum <- apply(abund_df, 1, sum)
      SmaxN_timestep <- max(vect_maxN_sum) 
      
      
      # return:
      return_list <- list("maxN" = maxN, 
                          "SmaxN_timestep" = SmaxN_timestep,
                          "maxN_cam" = maxN_cam,
                          "maxN_timestep" = maxN_timestep)
      
      return(return_list)
      
    } # end if fish speed is not taken into account

    
  } # if there is more than one camera
  
  # if there only one camera, return only maxN:
  
  if (ncol(dist_df) == 1) {
    print("SmaxN package has been computed for cases where more than one camera
          is used.")
    return(list("maxN"  = max(abund_df)))
  }  # if there is only one camera
  
  
}

