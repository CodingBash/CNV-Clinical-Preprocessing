
setwd("~/Documents/Git-Projects/Git-Research-Projects/drug-response-prediction")
source("helperFunctions.R")
source("cnProfileVisualizationLibrary.R")
source("facetsAnalysisLibrary.R")

setwd("~/Documents/Git-Projects/Git-Research-Projects/cn-ml-analysis")
source("featureMatrixAssignment.r")

setwd("~/Documents/Git-Projects/Git-Research-Projects/drug-response-prediction")
loaded_samples <- load_samples(classes = c("T", "F", "M"), sampleList = "sampleList.csv")

reference = "hN30"
dir = "output/FACETS_Reference_hN30_8_2_18_1/"
for(sample in loaded_samples){
    setwd("~/Documents/Git-Projects/Git-Research-Projects/FACETS_write_files")
    raw_segments <- retrieveFacetsSegments(sample, sample_subdir = "/", reference = reference, dir = dir)
    bed <- segmentsToBedFormat(raw_segments)

    setwd("~/Documents/Git-Projects/Git-Research-Projects/CNprep-Slicing-CORE-Analysis/")
    Acores <- retrieveCores("./output/coresResults/prev_run_8_2_2018_2/selectedCores/AselectedCoresBP.bed") # BED file of amplification recurrent regions
    Dcores <- retrieveCores("./output/coresResults/prev_run_8_2_2018_2/selectedCores/DselectedCoresBP.bed") # BED file of deletion recurrent regions
    ADcores <- retrieveCores("./output/coresResults/prev_run_8_2_2018_2/selectedCores/ADselectedCoresBP.bed") # BED file of both recurrent regions

    visualizeCNProfile(title = paste0("Profile for sample = ", sample), facets_segment_data = bed, Acores = Acores, Dcores = Dcores, save = FALSE)
}



