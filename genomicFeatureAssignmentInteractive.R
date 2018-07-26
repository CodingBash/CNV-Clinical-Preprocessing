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
cd_local("hT_results/finalCores/prev_run_7_3_2018")
Acores <- retrieveCores("AfinalCoresBP_subtractFraction.bed") # BED file of amplification recurrent regions
Dcores <- retrieveCores("DfinalCoresBP_subtractFraction.bed") # BED file of deletion recurrent regions

cd_local("resources")
aucData <- readRDS("listSampleTESAUC.RDS")


#
# Retrieve training set
#
cd_doc()
training_set <- retrieveTrainingSet(samples, Acores, Dcores, binDir = "CSHL/Project_TUV_12995_B01_SOM_Targeted.2018-03-02/")
training_set$matrix <- attachLabelsToSet(matrix_training_set = training_set$matrix, labelData = aucData)

visualizeUnclusteredHeatmap(training_set$melted)
hc <- clusterTrainingSet(training_set$melted, visualize = TRUE)
plot(hc)

cd_local("mlOutput")
write.csv(training_set$matrix, file ="coreTrainingSet_7_26_2018_1.csv")
