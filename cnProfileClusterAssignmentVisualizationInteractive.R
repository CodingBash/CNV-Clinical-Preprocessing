#
# Visualize CN profile that is color coded based on cluster assignment from CNprep
#

#
# This script produces the CN profile images and saves in preferred directory
#

#
# Load sources
#
setwd("~/Git-Projects/Git-Research-Projects/drug-response-prediction")
source("cnProfileClusterAssignmentVisualizationLibrary.R")

#
# Prepare parameters
#
cd_local()
loaded_samples <- load_samples(classes = c("T"), sampleList = "sampleList.csv")
chromosomeSizes <- readRDS("./chromosomeSizes.rds")
base_axis <- generateBaseAxis(chromosomeSizes)

for(sample in loaded_samples){
  visualizeClusterAssignmnet(sample, save = FALSE, saveDir = "./clusterAssignmentVisualization/")
}