#
# This script saves the GMM visualization for all hT samples
#

#
# Load library functions
# WARNING: If sources don't load, manually setwd() to src file location
#
source("helperFunctions.R")
source("visualizeGMMLibrary.R")

# Load names of hT samples
loaded_samples <- load_samples(classes = c("T"), sampleList = "sampleList.csv")

#
# Save GMM visualization
#
plots <- list()
for(i in seq(1, length(loaded_samples))){
  sample <- loaded_samples[[i]]
  
  #
  # Retrive CNprep::CNpreprocessing segtable output
  #
  cd_local()
  segtable <- retrieveSegtable(sample, dir = "segClusteringResults/prev_run3/")
  
  # Generate and save GMM visualization
  displayGMM(segtable = segtable, sample = sample, print = TRUE, save = FALSE)
}
