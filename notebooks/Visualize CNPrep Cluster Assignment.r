
setwd("~/Documents/Git-Projects/Git-Research-Projects/drug-response-prediction")
source("cnProfileClusterAssignmentVisualizationLibrary.R")

chromosomeSizes <- readRDS("./chromosomeSizes.rds")
base_axis <- generateBaseAxis(chromosomeSizes)

# Set the directory of the CNprep segtable results
segDir <- "segClusteringResults/prev_run1/"

# Set the sample that you want to view the CN profile for here.
loaded_samples <- c("hT1", "hT25","hT70", "hT72")

# With the 0-axis
for(sample in loaded_samples){
  visualizeClusterAssignment(sample, base_axis = base_axis, segDir = segDir, save = FALSE, saveDir = "./clusterAssignmentVisualization/")
}

# Without the 0-axis
for(sample in loaded_samples){
  visualizeClusterAssignment(sample, segDir = segDir, save = FALSE, saveDir = "./clusterAssignmentVisualization/")
}

loaded_samples <- load_samples(classes = c("T"), sampleList = "sampleList.csv")

for(sample in loaded_samples){
  visualizeClusterAssignment(sample, base_axis = base_axis, segDir = segDir, save = FALSE, saveDir = "./clusterAssignmentVisualization/")
}

for(sample in loaded_samples){
  visualizeClusterAssignment(sample, segDir = segDir, save = FALSE, saveDir = "./clusterAssignmentVisualization/")
}

visualizeClusterAssignment(sample = "hT1", categories = c("chr17"), base_axis = base_axis, segDir = segDir, save = FALSE, saveDir = "./clusterAssignmentVisualization/")


