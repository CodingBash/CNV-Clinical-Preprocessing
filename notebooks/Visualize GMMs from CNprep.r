
# WARNING: If sources don't load, manually setwd() to src file location
source("helperFunctions.R")
source("visualizeGMMLibrary.R")

sample_list <- c("hT3", "hT30") # Set the sample that you want to view the GMM for here.

for(sample in sample_list){
    cd_local()
    segtable_mclust_V <- retrieveSegtable(sample, dir = "segClusteringResults/prev_run3/")
    displayGMM(segtable = segtable_mclust_V, sample = sample, print = TRUE, save = FALSE)
}

for(sample in sample_list){
    cd_local()
    segtable_mclust_E <- retrieveSegtable(sample, dir = "segClusteringResults/prev_run1/")
    displayGMM(segtable = segtable_mclust_E, sample = sample, print = TRUE, save = FALSE)
}

