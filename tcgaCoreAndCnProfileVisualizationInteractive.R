#
# This script produces the CN profile images and saves in preferred directory
#

#
# Load sources
#
setwd("~/Git-Projects/Git-Research-Projects/drug-response-prediction")
source("helperFunctions.R")
source("cnProfileVisualizationLibrary.R")
source("facetsAnalysisLibrary.R")

tumorId = 6

cd_facets("output")
raw_segments <- retrieveFacetsSegmentsFromObject(tumorId, tumorId, fitPrefix = "facetsG5Fit_", dir = "prev_run_4/") # Get matching tumor-normal pair 
bed <- segmentsToBedFormat(raw_segments)

cd_facets("cores")
Acores <- read.table(paste0("AcoresBP_t", tumorId, ".bed"), header = FALSE, sep = "\t") 
Dcores <- read.table(paste0("DcoresBP_t", tumorId, ".bed"), header = FALSE, sep = "\t") 
visualizeCNProfile(facets_segment_data = bed, Acores = Acores, Dcores = Dcores, save = FALSE)

