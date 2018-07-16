
setwd("~/Documents/Git-Projects/Git-Research-Projects/drug-response-prediction")
source("helperFunctions.R")
source("visualizeGMMLibrary.R")

sample_list <- c("hT30") # Set the sample that you want to view the GMM for here.

for(sample in sample_list){
    setwd("~/Documents/Git-Projects/Git-Research-Projects/drug-response-prediction")
    segtable_mclust_V <- retrieveSegtable(sample, dir = "segClusteringResults/prev_run3/")
    displayGMM(segtable = segtable_mclust_V, sample = sample, print = TRUE, save = FALSE)
}

for(sample in sample_list){
    setwd("~/Documents/Git-Projects/Git-Research-Projects/drug-response-prediction")
    segtable_mclust_E <- retrieveSegtable(sample, dir = "segClusteringResults/prev_run1/")
    displayGMM(segtable = segtable_mclust_E, sample = sample, print = TRUE, save = FALSE)
}

for(sample in sample_list){
    setwd("~/Documents/Git-Projects/Git-Research-Projects/drug-response-prediction")
    segtable_mclust_V_minJoin_05 <- retrieveSegtable(sample, dir = "segClusteringResults/prev_run_7_15_2018_2/")
    displayGMM(segtable = segtable_mclust_E, sample = sample, print = TRUE, save = FALSE)
}

for(sample in sample_list){
    setwd("~/Documents/Git-Projects/Git-Research-Projects/drug-response-prediction")
    segtable_mclust_E_minJoin_05 <- retrieveSegtable(sample, dir = "segClusteringResults/prev_run_7_15_2018_3/")
    displayGMM(segtable = segtable_mclust_E, sample = sample, print = TRUE, save = FALSE)
}

for(sample in sample_list){
    setwd("~/Documents/Git-Projects/Git-Research-Projects/drug-response-prediction")
    segtable <- retrieveSegtable(sample, dir = "segClusteringResults/prev_run_7_15_2018_14/")
    displayGMM(segtable = segtable, sample = sample, print = TRUE, save = FALSE)
}

for(sample in sample_list){
    setwd("~/Documents/Git-Projects/Git-Research-Projects/drug-response-prediction")
    segtable <- retrieveSegtable(sample, dir = "segClusteringResults/prev_run_7_15_2018_15/")
    displayGMM(segtable = segtable, sample = sample, print = TRUE, save = FALSE)
}

for(sample in sample_list){
    setwd("~/Documents/Git-Projects/Git-Research-Projects/drug-response-prediction")
    segtable <- retrieveSegtable(sample, dir = "segClusteringResults/prev_run_7_15_2018_16/")
    displayGMM(segtable = segtable, sample = sample, print = TRUE, save = FALSE)
}

for(sample in sample_list){
    setwd("~/Documents/Git-Projects/Git-Research-Projects/drug-response-prediction")
    segtable <- retrieveSegtable(sample, dir = "segClusteringResults/prev_run_7_15_2018_17/")
    displayGMM(segtable = segtable, sample = sample, print = TRUE, save = FALSE)
}

for(sample in sample_list){
    setwd("~/Documents/Git-Projects/Git-Research-Projects/drug-response-prediction")
    segtable <- retrieveSegtable(sample, dir = "segClusteringResults/prev_run_7_15_2018_18/")
    displayGMM(segtable = segtable, sample = sample, print = TRUE, save = FALSE)
}

for(sample in sample_list){
    setwd("~/Documents/Git-Projects/Git-Research-Projects/drug-response-prediction")
    segtable <- retrieveSegtable(sample, dir = "segClusteringResults/prev_run_7_15_2018_19/")
    displayGMM(segtable = segtable, sample = sample, print = TRUE, save = FALSE)
}

for(sample in sample_list){
    setwd("~/Documents/Git-Projects/Git-Research-Projects/drug-response-prediction")
    segtable <- retrieveSegtable(sample, dir = "segClusteringResults/prev_run_7_15_2018_20/")
    displayGMM(segtable = segtable, sample = sample, print = TRUE, save = FALSE)
}


