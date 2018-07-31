#source("https://bioconductor.org/biocLite.R")
#biocLite("GenomicRanges")

setwd("~/Git-Projects/Git-Research-Projects/drug-response-prediction")
source("genomicFeatureAssignment.R")

#
# Load sample to retrieve feature set for
# TODO: There seems to be a scope conflict - samples is getting overwritten
#
samples <- load_samples(classes = c("T","F", "M"))

#
# Retrieve CORE features
#
setwd("~/Git-Projects/Git-Research-Projects/hN_core_artifacts")
ADcores <- retrieveCores("./hT_output/prev_run_7_30_2018_1/selectedCores/ADselectedCoresBP.bed") # BED file of amplification recurrent regions

cd_local("resources")
aucData <- readRDS("listSampleTESAUC.RDS")


#
# Retrieve training set
#
setwd("~/Git-Projects/Git-Research-Projects/FACETS_write_files")
training_set <- retrieveTrainingSet(loaded_samples = samples, ADcores = ADcores, sample_subdir = "/", reference = "hN31", dir = "output/FACETS_Reference_hN31_7_28_18_2/")
training_set$matrix <- attachLabelsToSet(matrix_training_set = training_set$matrix, labelData = aucData)

visualizeUnclusteredHeatmap(training_set$melted)
hc <- clusterTrainingSet(training_set$melted, visualize = TRUE)
plot(hc)

cd_local("mlOutput")
write.csv(training_set$matrix, file ="coreTrainingSet_7_31_2018_1.csv")
