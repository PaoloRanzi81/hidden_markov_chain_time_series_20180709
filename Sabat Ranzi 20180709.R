# # Improvements 20180706: 
# 1. bootstrapping; 
#   1.2. bootstrapping + MLE; 
#   1.3. It seems that the bootstrapped TMs decrease number of clusters as computation takes place. 
#        Possibly we should reset rEMM at each run; 
#        for resetting R at each bootstrapped run, try the following: rm(list = ls()) or 
#        closeAllConnections;
# 2. scaling: 
#   2.1. try scaled vs not-scaled version: does it make a difference?; 
#   2.2. robust scaling with R package WRS2; 
# 3. clustering results overlapped with ggplot2's results; 
# 5.thresholding: 
#   5.2. try different distances: eJaccard vs Euclidean;
#   5.3. try h_cluster, dendrogram and max function for selecting "a posteriori" the strongest
#        transition relationhips and strongest clusters; 
#   5.4. try whether removing and fading clusters could improve output; 
# 6. TFCE: since the output (i.e. transition probabilities matrix) is an adjacency matrix, 
#    TFCE could be used;
#   6.1. TFCE can handle only p-values or whatever estimator (like bootstrapping);  
# 7. centroids, medodis: how they work?
# 10. Implement PPV (e.g. false positives/false negatives ratio)
# 11. if necessary, implement test for Heterosckesasticity; 
# 12. if necessary, implement test for Multicollinearity; 
# 13. if necessary, implement test for Auto-correlation; 
# 14. if necessary, implement test for Stationarity; 
# 19. set.seed(33333)


###########################################################################################################
## DESCRIPTION: Hidden Markov Chain analyses of finacial time-series
## TUTORIAL for rEMM, see the following article:
# Hahsler, Michael, and Margaret H. Dunham. "remm: Extensible markov model for 
# data stream clustering in r." Journal of Statistical Software 35.5 (2010): 1-31.

## INSTRUCTIONS FOR THE R SCRIPT: 
# Caveats: 3.1. bootraapped is -unfortunately- only partially data-driven. Indeed, rEMM is dynamic: 
# the number of clusters could change with different runs of the bootstrapped dataset. 
# Therefore, select ONLY bootstrapped transition matrices from transition_matrix_resampled 
# with the number of clusters matching to the original. 



##########################################################################################################
## 1. LOAD PACKAGES
# 1.1. load libraries
library(tidyverse)
library(lubridate)
library(forecast)
library(rEMM)
library(Rgraphviz)
#library(boot)
#library(tictoc)


# 1.2. set pathway
set_directory<-("C:/Univ/scripts/R scripts/finance/sabat")
setwd(set_directory)


# 1.3. set  runs for bootstrapping: it is adviced to used 10000 runs;
#number_resamplings <- 10



##########################################################################################################
## 2. LOAD DATA
# 2.1. load file .csv 
# data_temp_1 <- read_csv("C:/time_series/sabat/binance--BAT-ETH.csv", col_names =F)
# data_temp_2 <- read_csv("C:/time_series/sabat/binance--EOS-BTC.csv", col_names =F)
data_temp_1 <- read_csv("C:/time_series/sabat/binance--BTC-USDT.csv", col_names =F)
data_temp_2 <- read_csv("C:/time_series/sabat/binance--LTC-USDT.csv", col_names =F)
# data_temp_3 <- read_csv("C:/time_series/sabat/binance--ETH-XLM.csv", col_names =F)


# 2.2. loop through either each day or each hour
for (aaa in c(1:10)){ 
  data_temp_sample <- data_temp_1 %>% 
  dplyr::select(c(X1, X4), -one_of("X3")) %>%
    # dplyr::select(-one_of("X3")) %>% 
  dplyr::filter(X1 >= (as.POSIXct("2018-06-29 22:00:03", tz = "UTC")) & X1 <= 
                  (as.POSIXct("2018-06-30 00:00:03", tz = "UTC"))) %>% as.data.frame() 
  # eliminate attributes from data.frame
  attr(data_temp_sample, c("spec", "row.names")) <- NULL
  
  data_temp_sample_1 <- data_temp_1 %>% 
    dplyr::select(c(X1, X4), -one_of("X3")) %>%
    # dplyr::select(-one_of("X3")) %>% 
    dplyr::filter(X1 >= (as.POSIXct("2018-06-29 22:00:03", tz = "UTC")) & X1 <= 
                    (as.POSIXct("2018-06-30 00:00:03", tz = "UTC"))) %>% as.data.frame() 
  # eliminate attributes from data.frame
  attr(data_temp_sample_1, c("spec", "row.names")) <- NULL

  data_temp_sample_2 <- data_temp_2 %>% 
    dplyr::select(c(X1, X4), -one_of("X3")) %>%
    # dplyr::select(-one_of("X3")) %>% 
    dplyr::filter(X1 >= (as.POSIXct("2018-06-29 22:00:03", tz = "UTC")) & X1 <= 
                    (as.POSIXct("2018-06-30 00:00:03", tz = "UTC"))) %>% as.data.frame() 
  # eliminate attributes from data.frame
  attr(data_temp_sample_1, c("spec", "row.names")) <- NULL
}


# 2.3. free up RAM memory
data_temp_sample <- NULL


# 2.4. delete uncessary objects from workspace
data_temp_1 <- NULL
data_temp_2 <- NULL
data_temp_3 <- NULL
# data_temp_4 <- NULL


# 2.5. concatenate time-series
data_temp <- data.frame(c(data_temp_sample_1, data_temp_sample_2[,2]))
# data_temp <- c(data_temp_sample_1, data_temp_sample_2[,2])
#data_temp <- full_join(data_temp_sample_1, data_temp_sample_2)
#data_temp <- left_join(data_temp_sample_1, data_temp_sample_2)


# 2.6. delete uncessary objects from workspace
data_temp_sample_1 <- NULL 
data_temp_sample_2 <- NULL 
data_temp_sample_3 <- NULL
# data_temp_sample_4 <- NULL


# 2.7. create dataset to be plugged in rEMM: such a step eliminate on purpose 
# the first column of the date
data_temp %>% select(-one_of(c("X1"))) -> data


## TEST: test scaling (e.g. by using log: it trasforms data set with huge variance 
# in one with less variance). In this case: it could be either center or not centered 
# scaling (in general: do not use scaling!)
# data_temp_scaled_centered <- scale(data_temp, center = TRUE)
# data_temp_scaled_not_centered <- scale(data_temp, center = FALSE)
# 
# data_temp <- data_temp_scaled_centered
# data_temp <- data_temp_scaled_not_centered 


## TEST: Linear regression 
# data_tempres <- lm(data_temp[[2]] ~  data_temp[[1]],  data=data_temp)
# summary(data_tempres)



##########################################################################################################
## 3. LOESS (PS: use "nrow" instead of first column as x-axis)
# 3.1. create time bins in the first column (useful for LOESS visualization)
data_temp$X1 <- c(1:nrow(data_temp))


# 3.2. compute LOESS
loess_credit <- loess(data_temp[[1]] ~ data_temp[[2]],  data=data_temp, span=0.05) # 10% smoothing span
smoothed_credit  <- predict(loess_credit)
loess_debt <- loess(data_temp[[1]] ~ data_temp[[2]],  data=data_temp, span=0.05) # 10% smoothing span
smoothed_debt <- predict(loess_debt)


# 3.3. plot LOESS
plot(y= data_temp[[2]], x=data_temp[[1]], type="p", col= "red", 
     main="Loess Smoothing and Prediction", xlab="date", ylab="USD min")
points(y= data_temp[[3]], x=data_temp[[1]], type="p", col= "blue" )
lines(smoothed_credit,  col="red")
lines(smoothed_debt,  col="blue")


# 3.4. zoom in time-series: the aim at sucha a step is to spot some otuliers by eyeballing the data
data_temp$X1 <- c(1:nrow(data_temp))
plot(y= data_temp[[2]], x=data_temp[[1]], type="l", col= "black", 
     main="binance--ETH-XLM.csv (subsample)", xlab="date", ylab="Bitcoin")
#plot(y= data_temp[[2]][1600:2002], x=data_temp[[1]][1600:2002], type="l", col= "black", 
#main="binance--EOS-BTC.csv (subsample)", xlab="date", ylab="USD min")
plot(y= data_temp[[3]][1600:2002], x=data_temp[[1]][1600:2002], type="l", col= "black", 
     main="~ 1 Hr: from 2018-05-15 03:21:37 to 2018-05-15 04:15:13", xlab="time", ylab="binance--EOS-BTC")
plot(y= data_temp[[4]][1600:2002], x=data_temp[[1]][1600:2002], type="l", col= "black", 
     main="Loess Smoothing and Prediction", xlab="date", ylab="USD min")


# TEST: load files again
#data_temp_1 <- read_csv("C:/time_series/sabat/binance--BAT-ETH.csv", col_names =F)
#data_temp_2 <- read_csv("C:/time_series/sabat/binance--EOS-BTC.csv", col_names =F)
#data_temp_3 <- read_csv("C:/time_series/sabat/binance--ETH-XLM.csv", col_names =F)



##########################################################################################################
## 4. ORIGINAL rEMM 
# 4.1. rule-of-thumb for computing a data-driven threshold (according to Michael Hahsler's tutorial). 
# At the same time trial-and-error is also possible in order to get the right amount of clusters. 
# The goal is to have from 4 to 8 clusters.
rows <- sample(1:nrow(data), 100)
emm_threshold_data_driven <- quantile(dist(data[rows,]), .25)
emm_threshold <- emm_threshold_data_driven
emm_threshold <- 0.0005


# 4.2. set type of clustering method (e.g. "eJaccard", "euclidean") 
emm <- EMM(threshold = emm_threshold, measure="eJaccard")


# 4.3. compute clusters
emm <- build(emm, data, verb=FALSE)


# 4.4. for diagnositcs: get clusters' details (by uncommenting rows below).  
cluster_counts(emm)


# 4.5. plot clusters. Methods could be: "graph". "MDS"; 
rEMM::plot(emm, method="graph")
rEMM::plot(emm, method="MDS")

# 4.6. double-check transition matrix (mainly for getting one-to-one score of probability)  
transition_matrix(emm, prior=FALSE) 
transition_matrix(emm, type="counts", prior=FALSE) 


# 4.7. remove self-transitions
emm_1 <- remove_selftransitions(emm, copy = TRUE)
emm<- emm_1
cluster_counts(emm_1)


# 4.8. are the following reliable centroids?
cluster_centers <- data.frame(cluster_centers(emm_1))


# TEST: dendrogram of clusters/centroids
# graphics::plot(x=ccc[,1], y=ccc[,2])
# graphics::plot(x=t(ccc[,1]), y=t(ccc[,2]))
# rEMM::plot(c, method="MDS")
# rEMM::plot(c, method="graph")
# fit <- hclust(c, method="ward") 
# data(ccc)
# class(cluster_centers)
# ddd<- transition_table(emm_1, data, type ="prob", initial_transition=TRUE) 
# eee<- data.frame(sort(ddd[,3]))


# 4.9. run again step 4.5. but now with updated data set
plot(emm_1, method="graph")
plot(emm_1, method="MDS")


# 4.10. run again step 4.6.
transition_matrix(emm_1, prior=FALSE) 
transition_matrix(emm_1, type="counts", prior=FALSE) 


# # 4.10. reclustering states (only to be used a few times, since the correct "threshold" 
# # should be found empirically i.e. by trial-and-error) [it does NOT work!]
# k <- c(2:6)
# emmc <- recluster_hclust(emm_1, k=k, h=NULL, prune=NULL, method ="average", copy = TRUE) 
# plot(attr(emmc, "cluster_info")$dendrogram)


# 4.11. pruning (namely, getting rid of small transition probabilities)
rare_threshold <- sum(cluster_counts(emm_1))*0.005
#plot(prune(emm_1, rare_threshold))
emm_2 <- prune(emm_1, rare_threshold)
transition_matrix(emm_1, prior=FALSE) 
#plot(emm_2, method="graph")


# 4.12. create transition matrix's object to be saved
transition_matrix_original <- transition_matrix(emm_2, prior=FALSE)


# 4.13. save matrix as .RData file
save(transition_matrix_original, file = file.path(set_directory, 
                                                  "transition_matrix_original.RData"))



##########################################################################################################
## 5. BOOTSTRAPPED rEMM (such a step is NOT mandatory! We can skip it and go directly to step 7 or step 8)
# 5.1. boostrapping (start of "for loop")
# 5.1.1. initialize transition_matrix_resampled (namely, all bootstrapped transition matrices)
transition_matrix_resampled <- NULL
# transition_matrix_output <- matrix(, nrow(data), ncol(data))

# TEST: start checking speed "for loop"
tic()


# 5.1.2. start "for loop" for bootstrapping
for (aaa in 1:number_resamplings){ 
indexes <- base::sample(c(1:nrow(data)), nrow(data), replace = TRUE)
data <- data[indexes,]


# 5.2. set details before computing clusters 
emm <- EMM(threshold=0.2, measure="eJaccard")
#emm <- EMM(threshold=1000, measure="euclidean")


# 5.3. compute clusters
emm <- build(emm, data,  verb=FALSE)


# 5.4. diagnositcs: get clusters' details (by uncommenting rows below), if extra pieces 
# of information are required 
cluster_counts(emm)


# 5.5. plot clusters
# plot(emm, method="graph")


# 5.6. double-check transition matrix (mainly for getting probability one-to-one scores)  
transition_matrix(emm, prior=FALSE) 
transition_matrix(emm, type="counts", prior=FALSE) 


# 5.7. remove self-transitions
emm_1 <- remove_selftransitions(emm, copy = TRUE)


# 5.8. run again 4.5. with updated data set
# plot(emm_1, method="graph")


# 5.9. run again 4.6.
transition_matrix(emm_1, prior=FALSE) 
transition_matrix(emm_1, type="counts", prior=FALSE) 


# # 5.10. reclustering states (only to be used a few times, since the correct "threshold"
# # should be found empirically i.e. by trial-and-error) [it does NOT work!]
# k <- c(2:6)
# emmc <- recluster_hclust(emm_1, k=k, h=NULL, prune=NULL, method ="average", copy = TRUE) 
# plot(attr(emmc, "cluster_info")$dendrogram)
# 
# sc <- sapply(emmc, score, EMMsim_test)
# names(sc) <- sq
# sc


# 5.11. pruning (=> getting rid of small transition probabilities)
rare_threshold <- sum(cluster_counts(emm_1))*0.005
rare_threshold
# plot(prune(emm_1, rare_threshold))
emm_2 <- prune(emm_1, rare_threshold)
# plot(emm_2, method="graph")


# 5.12. create transition matrix for each bootrapped dataset 
transition_matrix_like_original <- transition_matrix(emm_2, prior=FALSE)


# 5.13. Create "tidy table" (according to Headly Wickham's definition: the table 
# must be two dimensional and NOT multidimensional!) of transition_matrix_resampled: 
# all the single bootstrapped will be aggregated in one single R object
transition_matrix_resampled[aaa] <- list(transition_matrix_like_original)


# 5.14. boostrapping (end "for loop")
}


# TEST: end checking speed "for loop"
toc()


# 5.15. save transition_matrix_resampled as .RData file
save(transition_matrix_resampled, file = file.path(set_directory, 
                                                   "transition_matrix_resampled.RData"))



##########################################################################################################
# ## 6. QUANTILES
# # 6.1. rEMM is dynamic: the number of clusters could change with different runs of the 
# bootstrapped dataset. Therefore, select ONLY bootstrapped transition matrices from 
# transition_matrix_resampled with the number of clusters matching to the original. 
# transition_matrix_output_selected <- transition_matrix_output
# length(transition_matrix_output_selected) 
# for (bbb in 1:length(transition_matrix_resampled) {
#   transition_matrix_resampled_selected[bbb] <- length(
# as.vector(transition_matrix_resampled[[bbb]]) == length(as.vector(transition_matrix_original))
# }
# 
# 
# # 6.2. compute 95 % quintiles for each index of the selected bootstrapped transition matrices  
# for (ccc in 1:length(transition_matrix_output_selected) { 
#   for (ddd in 1:length(transition_matrix_output_selected)  # put the right indexing!
#   indexes <- base::sample(c(1:nrow(data)), nrow(data), replace = TRUE)
#   data <- data[indexes,]
#   
# 
# tot <- NULL
# for (eee in 1:length(transition_matrix_output)){ 
# tot[eee]<-transition_matrix_output[[g]][1,2]
# }
# quantile(tot, probs = c(0.025, 0.975), type = 8)
# 
# 
# }
# 
# class(df[,"Var1"])
# head(transition_matrix_output$Type)   



##########################################################################################################
# ## 7. THRESHOLDING OF transition_matrix_original
# #  The goal of such a step is building a comparison between transition_matrix_original 
# and the quantiles calculated from transition_matrix_resampled. Those quantiles should be 
# used for thresholding the transition_matrix_original, thus getting rid of the majoirity of 
# false positives. 
# #  7.1. load dataset transition_matrix_original
load(file.path(set_directory, "transition_matrix_original.RData"))



##########################################################################################################
## 8. PLOT: WARNING: SWITCH OFF THE FOLLOWING GRAPHICAL PART WHEN RUNNING BOOTSTRAPPING!!!
# 8.1. find visually clusters
catchment <- 2
plot(data[,catchment], type="l", ylab="scores",
     main=colnames(data)[catchment])


# 8.2. again find clusters
emm_2<-emm_1
state_sequence <- find_clusters(emm_2, data)


# 8.3. run again 4.3. with updated data set
#cluster_counts(emm_2)


# 8.4.Function for helping visualizing clusters
mark_states <- function(states, state_sequence, ys, col=0, label=NULL, ...) {
  x <- which(state_sequence %in% states)
  points(x, ys[x], col=col, ...)
  if(!is.null(label)) text(x, ys[x], label, pos=4, col=col)
}


# 8.5.Overlap raw data with the discovered clusters
mark_states("1", state_sequence, data[,catchment], col="blue", label="1")
mark_states("2", state_sequence, data[,catchment], col="red", label="2")
mark_states("3", state_sequence, data[,catchment], col="yellow", label="3")
mark_states("4", state_sequence, data[,catchment], col="purple", label="4")
mark_states("5", state_sequence, data[,catchment], col="green", label="5")
mark_states("6", state_sequence, data[,catchment], col="black", label="6")
mark_states("7", state_sequence, data[,catchment], col="violet", label="7")
mark_states("8", state_sequence, data[,catchment], col="orange", label="8")
mark_states("9", state_sequence, data[,catchment], col="brown", label="9")
mark_states("3", state_sequence, data[,catchment], col="purple", label="3")
mark_states("13", state_sequence, data[,catchment], col="blue", label="13")
mark_states("15", state_sequence, data[,catchment], col="red", label="15")
mark_states("16", state_sequence, data[,catchment], col="yellow", label="16")





