################################################################################
##
## Script to test the compute.SmaxN.bigUI() function which computes SmaxN for
## a given bloc (timestep + uniqueness interval)
##
## 1_test_computation_big_UI.R
##
## 27/04/2022
##
## Camille Magneville
##
################################################################################


# Exemple 1: Petit interval 
# https://docs.google.com/spreadsheets/d/1sNroVbL-gnZah0yV0ZbvYuM3jXbD1cioJKaQ9MCMEOM/edit#gid=0
# Table 1, time_df1

abund_df <- data.frame("A" = c(9,8,3,3,3,3,3), "B" = c(0,4,2,2,1,3,3), 
                       "C" = c(0,0,0,0,1,1,1), "D" = c(1,0,0,0,3,3,3))
  
time_df <- data.frame("A" = c(0,5,6,6), "B" = c(5,0,3,3), 
                      "C" = c(6,3,0,3), "D" = c(6,3,3,0))
rownames(time_df) <- c("A", "B", "C", "D")

timestep <- 3

value <- 4

SmaxN_small_UI <- 10 # computed before

compute.SmaxN.bigUI(abund_df, value, timestep, time_df,
                    SmaxN_small_UI) 


# Exemple 2: Exemple Seb: Big interval and have to return in cam
# Change timestep and value
# (https://docs.google.com/spreadsheets/d/1r613A9UOHpMxeg0Z02d-3yiM8Q_hxrtGFYvnXmLJZbY/edit#gid=0)
# Table 1, time_df1

abund_df <- data.frame("A" = c(9,8,3,3,3,3,3), "B" = c(0,4,2,2,1,3,3), 
                       "C" = c(0,0,0,0,1,1,1), "D" = c(1,0,0,0,3,3,3))

time_df <- data.frame("A" = c(0,4,5,5), "B" = c(4,0,2,2), 
                      "C" = c(5,2,0,2), "D" = c(5,2,2,0))
rownames(time_df) <- c("A", "B", "C", "D")

timestep <- 2

value <- 6

SmaxN_small_UI <- 10

compute.SmaxN.bigUI(abund_df, value, timestep, time_df,
                    SmaxN_small_UI) 




# Exemple 3: Exemple Seb: Big interval and have to return in cam
# Change timestep and value
# (https://docs.google.com/spreadsheets/d/1r613A9UOHpMxeg0Z02d-3yiM8Q_hxrtGFYvnXmLJZbY/edit#gid=0)

abund_df <- data.frame("A" = c(9,8,3,3,3,3,3), "B" = c(0,4,2,2,1,3,3), 
                       "C" = c(0,0,0,0,1,1,1), "D" = c(1,0,0,0,3,3,3))

time_df <- data.frame("A" = c(0,4,5,5), "B" = c(4,0,2,2), 
                      "C" = c(5,2,0,2), "D" = c(5,2,2,0))
rownames(time_df) <- c("A", "B", "C", "D")

timestep <- 1

value <- 6

SmaxN_small_UI <- 10

compute.SmaxN.bigUI(abund_df, value, timestep, time_df,
                    SmaxN_small_UI) 


# Exemple 4: Stopper des la premiere camera si sum trop petite
# https://docs.google.com/spreadsheets/d/1sNroVbL-gnZah0yV0ZbvYuM3jXbD1cioJKaQ9MCMEOM/edit#gid=0
# Table 2 Time_df3

abund_df <- data.frame("A" = c(10,0,0,0,0), "B" = c(0,0,0,0,1), 
                       "C" = c(0,0,0,0,3), "D" = c(0,0,0,0,1))

time_df <- data.frame("A" = c(0,1,3,3), "B" = c(1,0,2,2), 
                      "C" = c(3,3,0,2), "D" = c(3,2,2,0))
rownames(time_df) <- c("A", "B", "C", "D")

timestep <- 1

value <- 5

SmaxN_small_UI <- 10

compute.SmaxN.bigUI(abund_df, value, timestep, time_df,
                    SmaxN_small_UI) 


# Exemple 5: Stopper des la premiere camera si sum trop petite avec plusieurs ...
# ... valeurs a tester pour la premiere cam et plusieurs cell pour une valeur:
# https://docs.google.com/spreadsheets/d/1sNroVbL-gnZah0yV0ZbvYuM3jXbD1cioJKaQ9MCMEOM/edit#gid=0
# Table 3, Timedf3

abund_df <- data.frame("A" = c(10,0,0,0,0), "B" = c(0,0,1,1,4), 
                       "C" = c(0,0,0,0,3), "D" = c(0,0,0,0,1))

time_df <- data.frame("A" = c(0,1,3,3), "B" = c(1,0,2,2), 
                      "C" = c(3,2,0,2), "D" = c(3,2,2,0))
rownames(time_df) <- c("A", "B", "C", "D")

timestep <- 1

value <- 5

SmaxN_small_UI <- 10

compute.SmaxN.bigUI(abund_df, value, timestep, time_df,
                    SmaxN_small_UI)


# Exemple 6: question sur la premiere camera: change t on la premoere valeur

abund_df <- data.frame("A" = c(10,8,4,0), "B" = c(0,1,1,4), 
                       "C" = c(1,0,0,4), "D" = c(6,2,0,2))

time_df <- data.frame("A" = c(0,2,2,3), "B" = c(2,0,1,4), 
                      "C" = c(2,1,0,5), "D" = c(3,4,5,0))
rownames(time_df) <- c("A", "B", "C", "D")

timestep <- 1

value <- 4

SmaxN_small_UI <- 17

compute.SmaxN.bigUI(abund_df, value, timestep, time_df,
                    SmaxN_small_UI)
