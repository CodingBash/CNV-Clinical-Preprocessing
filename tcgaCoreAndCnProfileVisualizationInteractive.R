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




tumorId = 5

cd_facets("output")
raw_segments <- retrieveFacetsSegmentsFromObject(tumorId, tumorId, fitPrefix = "facetsG5Fit_", dir = "prev_run_4/") # Get matching tumor-normal pair 
bed <- segmentsToBedFormat(raw_segments)

cd_facets("cores")
Acores <- read.table(paste0("AcoresBP_t", tumorId, ".bed"), header = FALSE, sep = "\t") 
Dcores <- read.table(paste0("DcoresBP_t", tumorId, ".bed"), header = FALSE, sep = "\t") 
visualizeCNProfile(facets_segment_data = bed, Acores = Acores, Dcores = Dcores, categories = c("chr17"), save = FALSE)


#
# Iterate through each sample and generate CN profile
#
for(sample in loaded_samples){
    #
    # Retrieve and format input data
    #
    cd_doc()
    facets_segment_data <- retrieveFacetsSegments(sample, dir = "CSHL/Project_TUV_12995_B01_SOM_Targeted.2018-03-02/")
    facets_snp_data <- retrieveFacetsSnps(sample, dir = "CSHL/Project_TUV_12995_B01_SOM_Targeted.2018-03-02/")
    facets_segment_data <- segmentsToBedFormat(facets_segment_data)  
    facets_snp_data <- snpsToBedFormat(facets_snp_data)
    
    #
    # Generate CN profile from input data
    #
    cd_local()
    visualizeCNProfile(facets_segment_data = facets_segment_data, facets_snp_data = facets_snp_data, save = TRUE, saveDir = "CNprofiles/")
}
