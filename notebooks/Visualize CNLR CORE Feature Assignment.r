
setwd("~/Documents/Git-Projects/Git-Research-Projects/drug-response-prediction")
source("genomicFeatureAssignment.R")



#
# Load sample to retrieve feature set for
#
samples <- load_samples(classes = c("T", "M", "F"))

#
# Retrieve CORE features
#
setwd("~/Documents/Git-Projects/Git-Research-Projects/drug-response-prediction/hT_results/finalCores/prev_run_7_3_2018")
Acores <- retrieveCores("AfinalCoresBP_subtractFraction.bed") # BED file of amplification recurrent regions
Dcores <- retrieveCores("DfinalCoresBP_subtractFraction.bed") # BED file of deletion recurrent regions


options(warn=-1)

#
# Retrieve training set
#
setwd("~/Documents")
training_set <- retrieveTrainingSet(, Acores, Dcores, binDir = "CSHL/Project_TUV_12995_B01_SOM_Targeted.2018-03-02/")

visualizeUnclusteredHeatmap(training_set)

hc <- clusterTrainingSet(training_set, visualize = TRUE)
plot(hc)
