
setwd("~/Documents/Git-Projects/Git-Research-Projects/drug-response-prediction")
source("helperFunctions.R")
source("cnProfileVisualizationLibrary.R")
source("facetsAnalysisLibrary.R")

for(tumorId in seq(1, 7)){
    setwd(paste0("~/Documents/Git-Projects/Git-Research-Projects/FACETS_nonmatching_test/output"))
    raw_segments <- retrieveFacetsSegmentsFromObject(tumorId, tumorId, fitPrefix = "facetsG5Fit_", dir = "prev_run_4/") # Get matching tumor-normal pair 
    bed <- segmentsToBedFormat(raw_segments)

    setwd(paste0("~/Documents/Git-Projects/Git-Research-Projects/FACETS_nonmatching_test/cores"))
    Acores <- read.table(paste0("AcoresBP_t", tumorId, ".bed"), header = FALSE, sep = "\t") 
    Dcores <- read.table(paste0("DcoresBP_t", tumorId, ".bed"), header = FALSE, sep = "\t") 
    visualizeCNProfile(title = paste0("Profile for tumorId = ", tumorId), facets_segment_data = bed, Acores = Acores, Dcores = Dcores, save = FALSE)
}

