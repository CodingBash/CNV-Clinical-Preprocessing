
source("helperFunctions.R")
source("visualizeGMMLibrary.R")

sample <- "hT3" # Set the sample that you want to view the GMM for here.

cd_local()
segtable_mclust_V <- retrieveSegtable(sample, dir = "segClusteringResults/prev_run3/")
segtable_mclust_E <- retrieveSegtable(sample, dir = "segClusteringResults/prev_run1/")

displayGMM(segtable = segtable_mclust_V, sample = sample, print = TRUE, save = FALSE)

displayGMM(segtable = segtable_mclust_E, sample = sample, print = TRUE, save = FALSE)


