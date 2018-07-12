#
# This script contains a set of helper functions for the other CNA scripts
#

#
# Sets the directory to the project workspace
# TODO: directory is hardcoded
#
cd_local <- function() {
  setwd("C:/Users/bbece/Documents/Git-Projects/Git-Research-Projects/drug-response-prediction")
}

#
# Sets the directory to the project workspace
# TODO: directory is hardcoded
#
cd_doc <- function() {
  setwd("C:/Users/bbece/Documents")
}

#
# Sets the directory to HPC home directory
#
cd_home <- function() {
  setwd("~/code/hN_core_artifacts")
}
# 
# Load all sample names given a sample class
#
load_samples <- function(classes = c("N"), sampleList = "sampleList.csv") {
  samples <- read.table(sampleList, header=T, sep = "\t", stringsAsFactors = F)
  loaded_samples <- c(NA)
  loaded_samples.index <- 1
  for(sample in samples$Organoids){
    for(class in classes){
      if(substring(sample, 2,2) == class){
        loaded_samples[loaded_samples.index] <- sample
        loaded_samples.index <- loaded_samples.index + 1
        next
      }
    }
  }
  return(loaded_samples)
}