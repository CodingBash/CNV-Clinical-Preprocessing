
setwd("~/Documents/Git-Projects/Git-Research-Projects/drug-response-prediction")
source("helperFunctions.R")
source("visualizeGMMLibrary.R")

sample_list <- c("hT1", "hT25","hT70", "hT72") # Set the sample that you want to view the GMM for here.

for(sample in sample_list){
    setwd("~/Documents/Git-Projects/Git-Research-Projects/drug-response-prediction")
    segtable_mclust_E <- retrieveSegtable(sample, dir = "segClusteringResults/prev_run1/")
    displayGMM(segtable = segtable_mclust_E, sample = sample, print = TRUE, save = FALSE)
}

for(sample in sample_list){
    setwd("~/Documents/Git-Projects/Git-Research-Projects/drug-response-prediction")
    segtable_mclust_V <- retrieveSegtable(sample, dir = "segClusteringResults/prev_run3/")
    displayGMM(segtable = segtable_mclust_V, sample = sample, print = TRUE, save = FALSE)
}

for(sample in sample_list){
    setwd("~/Documents/Git-Projects/Git-Research-Projects/drug-response-prediction")
    segtable_mclust_V <- retrieveSegtable(sample, dir = "segClusteringResults/prev_run_7_19_2018_1/")
    displayGMM(segtable = segtable_mclust_V, sample = sample, print = TRUE, save = FALSE)
}

for(sample in sample_list){
    setwd("~/Documents/Git-Projects/Git-Research-Projects/drug-response-prediction")
    segtable_mclust_V <- retrieveSegtable(sample, dir = "segClusteringResults/prev_run_7_19_2018_2/")
    displayGMM(segtable = segtable_mclust_V, sample = sample, print = TRUE, save = FALSE)
}


