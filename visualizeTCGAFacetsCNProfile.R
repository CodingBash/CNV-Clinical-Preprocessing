setwd("~/Git-Projects/Git-Research-Projects/drug-response-prediction")
source("helperFunctions.R")
source("cnProfileVisualizationLibrary.R")
source("facetsAnalysisLibrary.R")

#
# Visualize all segments per each tumor overlapped among eachother
#
for(tumorId in seq(1, 7)){
  tumorId = 6 # TODO Hardcoded in
  all_tumor_segments <- data.frame()
  for(normalId in seq(1, 7)){
    #if(normalId == tumorId) next
    try({
    cd_facets()
    raw_segments <- retrieveFacetsSegmentsFromObject(tumorId, normalId, fitPrefix = "facetsG5Fit_", dir = "output/prev_run_1/")
    bed <- segmentsToBedFormat(raw_segments)
    all_tumor_segments <- rbind(all_tumor_segments, bed)
    }, silent = TRUE)
  }
  visualizeCNProfile(all_tumor_segments,  categories = c("chr17"), save = FALSE)
}





