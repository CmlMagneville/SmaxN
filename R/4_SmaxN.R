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
#' @param speed a numerical value refering to the maximal speed of the 
#'  studied species. \strong{Speed must be given in meters per second}. If the
#'  computation of maxN values should not take into account speed (that is
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



SmaxN.computation <- function(abund_df, speed, dist_df) {
  
  
  # correct name because fish_speed used in all the algo:
  fish_speed <- speed
  
  
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
  
 
  #### Begin the SmaxN search
  
  ### Checks if no easy solution (SmaxN_small_UI):
  
  # if there is more than  one camera:
  if (ncol(dist_df) != 1) {
    
    # if fish_speed is taken into account:
    if (! is.null(fish_speed)) {
      
      
      # get the time_df using the 1st function:
      time_df <- SmaxN::compute.cam.time(dist_df = dist_df, fish_speed = fish_speed)
      
      
      # get the length of the small interval build the lowest time between ...
      # ... camera pairs:
      small_UI <- min(apply(time_df[, ], 1, function(x) min(x[x > 0])))
      # idem with max:
      big_UI <- max(time_df)
      
      ## create one df for the false "SmaxN" values for each timestep for the ...
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
        small <- SmaxN::pseudoSmaxN.timestep(time_df = time_df,
                                        abund_df = abund_df,
                                        timestep = i,
                                        value = small_UI)
        
        # compute the SmaxN for the timestep and the big span:
        big <- SmaxN::pseudoSmaxN.timestep(time_df = time_df,
                                    abund_df = abund_df,
                                    timestep = i,
                                    value = big_UI)
        
        # add values in the small and big df:
        small_SmaxN_df <- dplyr::add_row(row = rownames(abund_df)[i], SmaxN = small, small_SmaxN_df)
        big_SmaxN_df <- dplyr::add_row(row = rownames(abund_df)[i], SmaxN = big, big_SmaxN_df)
        
        small_SmaxN_df$SmaxN <- as.numeric(small_SmaxN_df$SmaxN)
        big_SmaxN_df$SmaxN <- as.numeric(big_SmaxN_df$SmaxN)

      } # end computation of Smaxn for each timestep and small and big spans
      
      # remove first rows filled with Na and use numerical format:
      small_SmaxN_df <- small_SmaxN_df[-1, ]
      big_SmaxN_df <- big_SmaxN_df[-1, ]

      # Compute the max of SmaxN of small spans for all timesteps:
      max_small <- max(small_SmaxN_df$SmaxN)
      
      
      
      ## Remove timesteps not to study ie the one with big_UI SmaxN < max_small
      # and the associated rows in the mask:
      clean_big_SmaxN_df <- big_SmaxN_df[which(! big_SmaxN_df$SmaxN < max_small), ]
      
      # add the associated rows:
      add_clean_big_df <- clean_big_SmaxN_df
      row_base <- as.numeric(add_clean_big_df$row)
      row_to_add <- c()
      
      for (j in (1:length(row_base))) {
        
        # create a vector containing all possible rows to add ...
        # ... to reduce if needed thereafter if row_base + big_Ui > nrow(abund_df):
        possible_rows <- c(row_base[j]:(row_base[j] + big_UI))
        possible_rows <- possible_rows[which(possible_rows <= nrow(abund_df))]
        
        row_to_add <- append(row_to_add, possible_rows)
                             
      }
      
      # remove double rows:
      row_to_add <- unique(row_to_add)
      
      # complete the df: 
      for (j in row_to_add) {
        
        # if not already in the df:
        if (! j %in% row_base) {
          
          add_clean_big_df$row <- as.numeric(add_clean_big_df$row)
          add_clean_big_df <- dplyr::add_row(row = j, 
                                             SmaxN = as.numeric(big_SmaxN_df$SmaxN[which(big_SmaxN_df$row == j)]), 
                                             add_clean_big_df)
          
        }
      }
      
      # now add_clean_big_df contains the timesteps to study, in which the general SmaxN is
      
      
      # create a list which contains the list of timesteps not to study:
      # ... help to reduce the frame of possibles:
      trash_rows <- setdiff(rownames(abund_df), add_clean_big_df$row)
      
      
      ## Check is the SmaxN is already known (if for one timestep: ...
      # ... SmaxN big UI = SmaxN small UI = max (SmaxN big UI)):
      max_big <- max(big_SmaxN_df$SmaxN) 
      
      # span on the two dfs:
      for (k in small_SmaxN_df$row) {
        for (m in big_SmaxN_df$row) {
          
          # if the same row:
          if (k == m) {
            
            # if SmaxN is the same then stop, we have the SmaxN:
            if ((small_SmaxN_df$SmaxN[as.numeric(k)] == big_SmaxN_df$SmaxN[as.numeric(m)])) {
              if (small_SmaxN_df$SmaxN[as.numeric(k)] == max_big) {
                print("SmaxN small UI = SmaxN big UI")
                return(list("SmaxN" = small_SmaxN_df$SmaxN[which(small_SmaxN_df$row == k)],
                            "timestep" = k))
              }
            }
          }
        }
      }
      
      
      
      # order the rows by decreasing order so that study intervals with the ...
      # ... biggest SmaxN first:
      order_big_SmaxN_df <- dplyr::arrange(add_clean_big_df, dplyr::desc(SmaxN))
      order_SmaxN_df <- dplyr::arrange(order_big_SmaxN_df, dplyr::desc(SmaxN))
      
      
    #### order the abundance df so cameras ok:
      
      # compute the order in whch cameras should be viewed:
      cam_order <- SmaxN::cam.order(time_df)
      
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
        print(paste0(round((b/nrow(order_SmaxN_df))*100, 2), sep = "", "%"))
        
        
        
        # compute the frame of possible around the studied timestep:
        frame_possible_df <- SmaxN::frame.possible(T = as.numeric(order_SmaxN_df$row[b]), 
                                            time_df, abund_df, trash_rows)
        
        # reduce the frame_possible_df by deleting rows with all NAs:
        frame_possible_df <- frame_possible_df[rowSums(is.na(frame_possible_df)) != ncol(frame_possible_df), ]

        # remove columns with only 0 (but not first camera if the cell = 0):
        diff_0 <- colnames(frame_possible_df[, which(colSums(frame_possible_df, na.rm = TRUE) != 0)])
        
        
        # if only one camera with values different from 0, then SmaxN equal the max of the cam:
        if (is.null(diff_0)) {
          
          # add the name of the camera for with the SmaxN is obtained (because only one cam with values != 0)
          if (v == max_small) {
            # name of the cam instead of the path:
            path_saved <- append(path_saved, names(which(colSums(frame_possible_df, na.rm = TRUE) != 0) == TRUE))
            # value:
            v <- max(frame_possible_df[, which(colSums(frame_possible_df, na.rm = TRUE) != 0)])
            # add the SmaxN of the given timestep to the SmaxN vect:
            SmaxN_vect <- append(SmaxN_vect, v)
          }
          
          if (v > max_small) {
            # value:
            v <- max(frame_possible_df[, which(colSums(frame_possible_df, na.rm = TRUE) != 0)])
            # add the SmaxN of the given timestep to the SmaxN vect:
            SmaxN_vect <- append(SmaxN_vect, v)
            # replace SmaxN:
            max_small <- v
            # name of the cam instead of the path:
            path_saved <- NULL
            path_saved <- names(which(colSums(frame_possible_df, na.rm = TRUE) != 0) == TRUE)
          }
          
          # do nothing if v < max_small because not interesting
          
          b <- b + 1
          
        } # end if only one camera with values different from 0
        
        
        # if more than one camera has values different from 0:
        if (! is.null(diff_0)) {
          
          # remove columns with all 0:
          if (colnames(frame_possible_df)[1] %in% diff_0) {
            frame_possible_df <- frame_possible_df[, diff_0]
          }
          if (! colnames(frame_possible_df)[1] %in% diff_0) {
            frame_possible_df <- frame_possible_df[, c(colnames(frame_possible_df)[1], diff_0)]
          }
          
          # if frame of possible is not filled with 0 (ie if some camera are different from all 0):
          if (is.data.frame(frame_possible_df)) {
            
            # delete rows with NA again (because some cameras may have been deleted):
            frame_possible_df <- frame_possible_df[rowSums(is.na(frame_possible_df)) != ncol(frame_possible_df), ]
            
            
            # check that with the first cell of this timestep we can have a "good"
            # SmaxN:
            
            first_cam_nm <- colnames(frame_possible_df)[1]
            
            first_cell_df <- data.frame("values" = abund_df[as.numeric(order_SmaxN_df$row[b]), 
                                                            first_cam_nm],
                                        "cam_nm" = first_cam_nm,
                                        "timestep" = as.numeric(order_SmaxN_df$row[b]))
            
            
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
              path_df <-  as.data.frame(matrix(ncol = 3, nrow = ncol(frame_possible_df)))
              colnames(path_df) <- c("value", "cam_nm", "timestep")
              
              # add the first cam:
              path_df[1, ] <- c(first_cell_df2$values, 
                                first_cell_df2$cam_nm,
                                first_cell_df2$timestep)
              
              
              #### First path to work on : 
              # the one on the studied timestep:
              
              for (j in (2:nrow(path_df))) {
                
                cam_nm <- colnames(frame_possible_df)[j]
                timestep <- first_cell_df2$timestep
                value <- frame_possible_df[which(rownames(frame_possible_df) == timestep), cam_nm]
                
                # complete the path_df
                path_df[j, ] <- c(value, cam_nm, timestep)
                
              }
              
              # Compute the sum for this first path and compare to SmaxN_small_UI:
              path_df$value <- as.numeric(path_df$value)
              S <- sum(path_df$value, na.rm = TRUE)
              
              if (S > max_small) {
                max_small <- S
              }
              
              # remove the last camera of the path,thus initialising the
              # recursive fct and running backward:
              path_df[nrow(path_df), ] <- rep(NA, 3)
              
              #### Compute all possible paths and keep in mind the highest SmaxN:
              
              # compute the SmaxN of the big interval of the given timestep:
              T <- rownames(frame_possible_df[!is.na(frame_possible_df[, 1]), ])
              T <- as.numeric(T)
              
              #print("Starting recursive ")
              
              value <- SmaxN::recursive.paths(T = T, 
                                       frame_possible_df, 
                                       n = ncol(frame_possible_df), 
                                       path_df = path_df, 
                                       SmaxN_small_UI = max_small,
                                       time_df)
              v <- value$"SmaxN_timestep"
              
              #print("Ending recursive ")
              
              # add the SmaxN of the given timestep to the SmaxN vect:
              SmaxN_vect <- append(SmaxN_vect, v)
              
              # add the path:
              if (v == max_small) {
                path_saved <- append(path_saved, list(value$"SmaxN_path"))
              }
              
              if (v > max_small) {
                max_small <- v
                path_saved <- NULL
                path_saved <- c(list(value$"SmaxN_path"))
              }
              
              b <- b + 1
              
            } # end if first cell possibly leads to a "good" SmaxN
            
          } # end if there is some colonnes different 0 in the frame of possible
          
          # if all columns in the frame of possible = 0:
          if (is.numeric(frame_possible_df)) {
            b <- b + 1
          }
          
        } # end if more than one camera have values different from 0
        

      } # end while
      
      
      # Compute the general SmaxN:
      SmaxN <- max(SmaxN_vect)
      
      # Compute other metrics:
      maxN <- max(abund_df)
      maxN_cam <- apply(abund_df, 2, max)
      maxN_timestep <- apply(abund_df, 1, max)
      
      vect_maxN_sum <- apply(abund_df, 1, sum)
      SmaxN_timestep <- max(vect_maxN_sum) 
      
      # removing NULL element from path_saved (the first one when creating object):
      path_saved <- path_saved[which(! is.null(path_saved))]
      
      
      # return:
      return_list <- list("maxN" = maxN, 
                          "SmaxN" = SmaxN,
                          "path_saved" = path_saved,
                          "number_SmaxN_path" =  length(path_saved),
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

