#
# This script produces the CN profile images and saves in preferred directory
#

#
# Load sources
# WARNING: If sources don't load, manually setwd() to src file location
#
setwd("C:/Users/bbece/Documents/Git-Projects/Git-Research-Projects/drug-response-prediction")
source("cnProfileVisualizationLibrary.R")
source("helperFunctions.R")

# Load samples to generate CN profiles for
loaded_samples <- load_samples(classes = c("N", "T"), sampleList = "sampleList.csv")

#
# Iterate through each sample and generate CN profile
#
for(sample in loaded_samples){
    #
    # Retrieve and format input data
    #
    cd_doc()    
    facets_snp_data <- retrieveFacetsSnps(sample, dir = "CSHL/Project_TUV_12995_B01_SOM_Targeted.2018-03-02/")
    facets_snp_data <- snpsToBedFormat(facets_snp_data)
    
    cd_local()
    facets_segment_clusters <- retrieveFacetsSegmentClusters(sample, dir = "segClusteringResults/prev_run1/")
    facets_segment_clusters <- segmentClustersToBedFormat(facets_segment_clusters)
      
    #
    # Generate CN profile from input data
    #
    cd_local()
    visualizeCNProfile(facets_segment_data = facets_segment_data, facets_snp_data = facets_snp_data, save = TRUE, saveDir = "CNprofilesUsingClusters/")
}
