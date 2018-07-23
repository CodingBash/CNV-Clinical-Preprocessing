
setwd("~/Documents/Git-Projects/Git-Research-Projects/drug-response-prediction")
source("cnProfileVisualizationLibrary.R")
source("helperFunctions.R")


segDir <- "segClusteringResults/prev_run1/"

# Set the sample that you want to view the CN profile for here.
loaded_samples <- c("hT1", "hT25","hT70", "hT72")

#
# Iterate through each sample and generate CN profile
#
for(sample in loaded_samples){
    setwd("~/Documents/Git-Projects/Git-Research-Projects/drug-response-prediction")
    facets_segment_clusters <- retrieveFacetsSegmentClusters(sample, dir = segDir)
    facets_segment_clusters <- segmentClustersToBedFormat(facets_segment_clusters)
      
    #
    # Generate CN profile from input data
    #
    visualizeCNProfile(facets_segment_data = facets_segment_clusters, facets_snp_data = facets_snp_data, save = FALSE)
}


# Load samples to generate CN profiles for
loaded_samples <- load_samples(classes = c("T"), sampleList = "sampleList.csv")

#
# Iterate through each sample and generate CN profile
#
for(sample in loaded_samples){
    #
    # Retrieve and format input data
    #
    setwd("~/Documents/Git-Projects/Git-Research-Projects/drug-response-prediction")
    facets_segment_clusters <- retrieveFacetsSegmentClusters(sample, dir = segDir)
    facets_segment_clusters <- segmentClustersToBedFormat(facets_segment_clusters)
      
    #
    # Generate CN profile from input data
    #
    cd_local()
    visualizeCNProfile(facets_segment_data = facets_segment_clusters, facets_snp_data = facets_snp_data, save = FALSE)
}

