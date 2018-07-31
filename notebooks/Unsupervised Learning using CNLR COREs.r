
#source("https://bioconductor.org/biocLite.R")
#biocLite("GenomicRanges")
setwd("~/Documents/Git-Projects/Git-Research-Projects/drug-response-prediction")
source("genomicFeatureAssignment.R")

#
# Load sample to retrieve feature set for
# TODO: There seems to be a scope conflict - samples is getting overwritten
#
setwd("~/Documents/Git-Projects/Git-Research-Projects/drug-response-prediction")
samples <- load_samples(classes = c("T","F", "M"), sampleList = "sampleList.csv")

#
# Retrieve CORE features
#
setwd("~/Documents/Git-Projects/Git-Research-Projects/hN_core_artifacts")
ADcores <- retrieveCores("./hT_output/prev_run_7_30_2018_1/selectedCores/ADselectedCoresBP.bed") # BED file of amplification recurrent regions

head(ADcores)

setwd("~/Documents/Git-Projects/Git-Research-Projects/drug-response-prediction")
aucData <- readRDS("./resources/listSampleTESAUC.RDS")

head(aucData$Gemcitabine)
head(aucData$Paclitaxel)
head(aucData$`SN-38`)
head(aucData$`5-FU`)
head(aucData$Oxaliplatin)

#
# Retrieve training set
#
setwd("~/Documents/Git-Projects/Git-Research-Projects/FACETS_write_files")
training_set <- retrieveTrainingSet(loaded_samples = samples, ADcores = ADcores, sample_subdir = "/", reference = "hN31", dir = "output/FACETS_Reference_hN31_7_28_18_2/")
training_set$matrix <- attachLabelsToSet(matrix_training_set = training_set$matrix, labelData = aucData)

options(repr.plot.width=15, repr.plot.height=15)
visualizeUnclusteredHeatmap(training_set$melted)

options(repr.plot.width=15, repr.plot.height=15)
hc <- clusterTrainingSet(training_set$melted, visualize = TRUE)

options(repr.plot.width=15, repr.plot.height=7)
plot(hc)

options(repr.plot.width=15, repr.plot.height=15)


for(label.i in seq(6)){
    labelData <- aucData[[label.i]]
    test_set <- cbind(training_set$melted)
    
    test_set$sampleId <- sapply(test_set$sampleId, 
                                function(sample){
                                    return(paste0(sample, "-", label.i, "-", labelData[labelData$SampleId == sample, ]$AUC))
                                })
    clusterTrainingSet(test_set, visualize = TRUE)
}



#setwd("~/Documents/Git-Projects/Git-Research-Projects/drug-response-prediction/mlOutput")
#write.csv(training_set$matrix, file ="coreTrainingSet_7_31_2018_1.csv")
