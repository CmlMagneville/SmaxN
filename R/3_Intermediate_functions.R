################################################################################
##
## 3_Intermediate_functions.R
##
## This script gathers three functions used in the global SmaxN function 
## (script 5): 
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
#' @noRd
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
#' @noRd
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
      # values <- c()
      # 
      # # loop on camera already ok:
      # for (j in (1:length(cam_ok))) {
      #   cam_done <- cam_ok[j]
      #   span <- time_df[cam_done, colnames(abund_df)[i]]
      #   values <- append(values, span)
      # } # end loop to get distances btw studied cam i and cam already ok
      # 
      # # get the span value which is:
      # # min( min(dist_cam_i_all_cam_ok) , min(dist_cam_i_cam1))
      # # = min(values):
      # span <- min(values)
      # span_start <- T - span
      # span_end <- T + span
      
      # the span is equal to the distance to the first camera:
      span <- time_df[cam1, colnames(abund_df)[i]]
      span_start <- T - span
      span_end <- T + span
      
      cam_ok <- append(cam_ok, colnames(abund_df)[i])
      
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
    if (i > 1) {
      
      if (span_start == 1 & span_end < nrow(abund_df)) {
        rownames(abund_df) <- as.numeric(rownames(abund_df))
        abund_df[which(as.numeric(rownames(abund_df)) > span_end), i] <- NA
      }
      
      if (span_start > 1 & span_end == nrow(abund_df)) {
        rownames(abund_df) <- as.numeric(rownames(abund_df))
        abund_df[which(as.numeric(rownames(abund_df)) < span_start), i] <- NA
      }
      
      if (span_start > 1 & span_end < nrow(abund_df)) {
        rownames(abund_df) <- as.numeric(rownames(abund_df))
        abund_df[which(as.numeric(rownames(abund_df)) > span_end), i] <- NA
        abund_df[which(as.numeric(rownames(abund_df)) < span_start), i] <- NA
      }
      
    }
    
    
    
  } # end loop on all cameras
  
  return("frame_possible" = abund_df)
  
}







#' Compute the biggest SmaxN for all the possible path of a given timestep
#'
#' This function computes all the paths beggining with a given timestep of the
#' first camera and return the highest SmaxN found on those path.
#' 
#' 
#' @param T a numerical value referring to the number of the row of the 
#' abundance dataframe studied (timestep)
#' 
#' @param frame_possible_df a dataframe which represents the frame of possible for
#' the studied timestep. Computed using the \code{frame.possible} function.
#' 
#' @param n a numerical value refering to the studied camera on which the
#' function searches new values
#' 
#' @param path_df a dataframe with value/camera_name/timestep information of 
#' cells already on the path
#' 
#' @param SmaxN_small_UI a numeric value refering to the highest pseudo SmaxN
#' found on the lowest intervals.
#'
#' @param time_df a numerical dataframe containing the minimal time 
#'  needed for an individual of the studied species to go from a camera to 
#'  another camera.There are as many rows as there are cameras and there are
#'  as many columns as there are cameras, thus the dataframe is symmetrical
#'  and the diagonal is filled with 0. This dataframe is the output of the 
#'  \code{compute.cam.time} function
#' 
#' @return the highest SmaxN value found on all the possble path given the 
#' studied timestep and the frame of possible.
#'
#' @noRd
#'  


recursive.paths <- function(T, frame_possible_df, n, path_df, SmaxN_small_UI, time_df) {
  
  
  #print(paste0("n = ", sep = "", n))
  
  
  # save the frame_possible df so can rewrite values when set to NA:
  saved_frame_possible_df <- frame_possible_df
  
  
  # get the partial sum pS which is the sum of the values alreday on the path:
  path_df$value <- as.numeric(path_df$value)
  pS <- sum(path_df$value, na.rm = TRUE)
  
  # initialise the vector saving the path for the highest SmaxN:
  path_saved <- NULL
  
  # while all cameras not studied:
  while (n > 1) {
    
    #print(paste0("n = ", sep = "", n))
    
    #### STILL SOME POSSIBLE CELLS FOR CAM N
    
    # while there are still some possible cells in the frame of possible for the camera n:
    while(length(frame_possible_df[which(! is.na(frame_possible_df[, n])), n]) > 0) {
      
      
      #print(paste0("Still some possible cells for cam", sep = " ", n))
      
      # if the camera is the last one, chose the maximal possible value:
      if (n == ncol(frame_possible_df)) {
        
        # check that the path so far is ok, if not we can just check the (n-1) cam:
        S <- sum(path_df$value, na.rm = TRUE) + 
          max(frame_possible_df[, n], na.rm = TRUE)
        
        if (S <= SmaxN_small_UI) {
          #print("raccourci")
          
          # this cell on (n-1) cam become already tested in frame_of_possible df:
          cam <- path_df$cam_nm[n - 1]
          timestep <- as.numeric(path_df$timestep[n - 1])
          frame_possible_df[which(rownames(frame_possible_df) == as.character(timestep)), cam] <- NA
          
          # the cell chosen for the (n-1) cam does not lead to good SmaxN:
          path_df[(n - 1), ] <- rep(NA, 3)
          
          n <- n - 1
          break
        }
        
        if (S > SmaxN_small_UI) {
        
          # check which cells are really possible given the path:
          cells_coord <- as.data.frame(matrix(ncol = 3, 
                                              nrow = length(frame_possible_df[which(! is.na(frame_possible_df[, n])), n])))
          colnames(cells_coord) <- c("values", "cam_nm", "timestep")
          cells_coord$cam_nm <- rep(colnames(frame_possible_df)[n], nrow(cells_coord))
          cells_coord$timestep <- rownames(frame_possible_df[which(! is.na(frame_possible_df[, n])), ])
          cells_coord$values <- frame_possible_df[which(! is.na(frame_possible_df[, n])), n]
          
          
          # check if cell is really possible given the path so far:
          possible_cell <- next.possible(max_bloc = path_df, 
                                         next_possible = cells_coord, 
                                         time_df)
          
          # reduce the cells_coord to the possible only (there is at least one possible cell):
          cells_coord_poss <- possible_cell[which(possible_cell$possible == TRUE), ]
          cells_coord_poss <- cells_coord_poss[, -ncol(cells_coord_poss)]
          
          #print(paste0("Cells_coord_poss = ", sep = "", cells_coord_poss))
          
          pS <- sum(path_df$value, na.rm = TRUE)
          value <- max(cells_coord_poss$values)
          S <- pS + value
          #print(paste0("S value = ", sep ="", S))
          #print(paste0("pS = ", sep ="", pS))
          #print(paste0("value = ", sep ="", value))
          
          
          
          # for the while loop to end, add NA in every cell of the frame of possible
          # ... for the studied cam:
          frame_possible_df[which(! is.na(frame_possible_df[, n])), n] <- NA
          
          #print(paste0("path so far is = ", sep = "", path_df))
          
          # add the value to the path:
          path_df$value[n] <- value
          path_df$cam_nm[n] <- colnames(frame_possible_df)[n]
          interm_timestep <- cells_coord_poss[which(cells_coord_poss$values == value), ]
          path_df$timestep[n] <- interm_timestep$timestep[1]
          path_df$timestep <- as.numeric(path_df$timestep)
          
          
          # if S > SmaxN_small_UI, replace (= so have at least one path on path_saved):
          # (= so have at least one path on path_saved):
          if (S == SmaxN_small_UI) {
            path_saved <- append(path_saved, list(path_df))
          }
          if (S > SmaxN_small_UI) {
            SmaxN_small_UI <- S
            path_saved <- 0
            path_saved <- list(path_df)
          }
          
          #print(paste0("path so far is = ", sep = "", path_df))
          #print(paste0("SmaxN_small_UI =", sep = "", SmaxN_small_UI))
          
        } # end if path interesting
           
      } # end if the camera is the last one
      
      
      # if the camera is not the last one, check all paths:
      if (n < ncol(frame_possible_df)) {
        
        # check that the path so far is ok, if not we can just check the (n-1) cam:
        S <- sum(path_df$value, na.rm = TRUE) + 
          sum(apply(frame_possible_df[, c(n:ncol(frame_possible_df))], 2, max, na.rm = TRUE))
        
        if (S <= SmaxN_small_UI) {
          
          if (n %in% c(2, 3, 4, 5)) {
            #print(paste0("n =", sep = " ", n))
            #print(paste0("timestep (n - 1) is", sep = " ", path_df$timestep[n-1]))
            #print(paste("SmaxN_small = ", SmaxN_small_UI))
            #print("S <= SmaxN_small")
          }
          
          # this cell on (n-1) cam become already tested in frame_of_possible df:
          cam <- path_df$cam_nm[n - 1]
          timestep <- as.numeric(path_df$timestep[n - 1])
          frame_possible_df[which(rownames(frame_possible_df) == timestep), cam] <- NA
          
          # give back the values to the cell of the (n) cam in case there is some NA
          # ... due to the "raccourci":
          # as we just want not to test paths with the cell at the (n-1) cam:
          frame_possible_df[, n] <- saved_frame_possible_df[, n]
          
          # the cell chosen for the (n-1) cam does not lead to good SmaxN:
          path_df[(n - 1), ] <- rep(NA, 3)
          
          n <- n - 1
          
          break
        }
        
        if (S > SmaxN_small_UI) {
          
          if (n %in% c(2, 3, 4, 5)) {
            #print(paste0("n =", sep = " ", n))
            #print(paste0("timestep (n - 1) is", sep = " ", path_df$timestep[n-1]))
            # print(paste("SmaxN_small = ", SmaxN_small_UI))
            # print("S > SmaxN_small")
          }
        
          # get coord of the first cell in the frame of possible:
          cells_coord <- as.data.frame(matrix(ncol = 3, 
                                              nrow = 1))
          colnames(cells_coord) <- c("values", "cam_nm", "timestep")
          cells_coord$cam_nm <- colnames(frame_possible_df)[n]
          cells_coord$timestep <- rownames(frame_possible_df[which(! is.na(frame_possible_df[, n])), ])[1]
          cells_coord$values <- frame_possible_df[which(! is.na(frame_possible_df[, n])), n][1]
          
          
          # check if cell is really possible given the path so far:
          possible_cell <- next.possible(max_bloc = path_df, 
                                         next_possible = cells_coord, 
                                         time_df)
          
          # reduce the cells_coord to the possible only (there is at least one possible cell):
          cells_coord_poss <- possible_cell[which(possible_cell$possible == TRUE), ]
          cells_coord_poss <- cells_coord_poss[, -ncol(cells_coord_poss)]
          
          #print(paste0("The studied cell is =", sep = "", possible_cell))
          
          
          # if the cell is not possible given the path:
          if (possible_cell$possible == FALSE) {
            
            # the cell becomes NA:
            cam_nm <- possible_cell$cam_nm
            timestep <- possible_cell$timestep
            frame_possible_df[which(rownames(frame_possible_df) == timestep), cam_nm] <- NA
            
            #print(paste0("FALSE - The frame of possible so far is =", sep = "", frame_possible_df))
            
          } # end if the cell is not possible given the path
          
          
          # if the cell is possible given the path:
          if (possible_cell$possible == TRUE) {
            
            # Compute the sum to compare with SmaxN_small_UI:
            path_df$value <- as.numeric(path_df$value)
            pS <- sum(path_df$value, na.rm = TRUE)
            
            # if the second to last camera (because only one camera left to add)
            if (n == (ncol(frame_possible_df) - 1)) {
              S <- pS + possible_cell$values + max(frame_possible_df[, c((n+1):ncol(frame_possible_df))], na.rm = TRUE)
            }
            
            # if any other camera:
            if (n < (ncol(frame_possible_df) - 1)) {
              S <- pS + possible_cell$values + 
                sum(apply(frame_possible_df[, c((n+1):ncol(frame_possible_df))], 2, max, na.rm = TRUE))
            }
            
            # if S > SmaxN_small_UI, ok can take the cell in the path (don't change because don't know if path possible):
            if (S > SmaxN_small_UI) {
              
              #print("S > SmaxN_small_UI")
              
              # add the cell to the path:
              path_df$value[n] <-  possible_cell$values
              path_df$cam_nm[n] <- cells_coord_poss$cam_nm
              path_df$timestep[n] <- cells_coord_poss$timestep
              path_df$timestep <- as.numeric(path_df$timestep)
              
              #print("Begin the recursive fct for n + 1 cam")
              #print(paste0("The frame of possible before is =", sep = "", frame_possible_df))
              
              # run for (n+1):
              n <- n + 1
              
              # this cell on n become already tested:
              frame_possible_df[cells_coord_poss$timestep, cells_coord_poss$cam_nm] <- NA
              
              #print(paste0("The frame of possible after is =", sep = "", frame_possible_df))
              #print(paste0("The path so far is ", sep = "", path_df))
              
              break
              
            } # end if S > SmaxN_small_UI
            
            
            # if S <= SmaxN_small_UI, abort mission: cell cam n becomes NA 
            # ... and give back the values of cam (n+1):
            if (S <= SmaxN_small_UI) {
              
              #print("S <= SmaxN_small_UI")
              
              # this cell on n become already tested:
              frame_possible_df[which(rownames(frame_possible_df) == possible_cell$timestep), possible_cell$cam_nm] <- NA
              
              # give back the values for cam (n+1) until nb_cam (n+1, n+2, n+3 ...):
              for (k in ((n+1):ncol(frame_possible_df))) {
                frame_possible_df[, k] <- saved_frame_possible_df[, k]
              }
              
              
            } # end if S > SmaxN_small_UI
          }
        
        } # end if the path is interesting
        
      } # end if not the last camera
      
      
    } # end while there are still some possible cells for the camera n
    
    
    
    #### NO MORE POSSIBLE CELLS FOR CAM N:
    
    if (length(frame_possible_df[which(! is.na(frame_possible_df[, n])), n]) == 0) {
      
      #print(paste0("No more possible cells for cam", sep = " ", n))
      
      # if we study the camera 2 (the timestep has been finished to be studied):
      if (n == 2) {
        
        n <- n - 1
        
      } # end if we study cam2
      
      
      # if we don't study cam2:
      if (n > 2) {
        
        # the cell from cam (n-1) becomes NA
        cam_nm <- path_df$cam_nm[(n-1)]
        timestep <- path_df$timestep[(n-1)]
        frame_possible_df[which(rownames(frame_possible_df) == timestep), cam_nm] <- NA
        
        #print(paste0("The frame possible df = ", sep = " ", frame_possible_df))
        
        # give back the values for the cells cam n:
        frame_possible_df[, n] <- saved_frame_possible_df[, n]
        
        #print(paste0("The frame possible df = ", sep = " ", frame_possible_df))
        
        # update path and remove n and n - 1 cam:
        path_df[c(n, n-1), ] <- rep(NA, 6)
        
        #print(paste0("The path df = ", sep = " ", path_df))
        
        # and then get interested about n - 1 cam:
        n <- n - 1
        
      } # end if we don't study cam2
      
      
    } # end if no more possible cells for cam n
    
    
  } # end while n > 1
  
  
  return(list("SmaxN_timestep"= SmaxN_small_UI,
              "SmaxN_path" = path_saved))
  
}
