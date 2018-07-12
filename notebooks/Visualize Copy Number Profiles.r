
source("https://bioconductor.org/biocLite.R")
biocLite("gtrellis")
biocLite("ComplexHeatmap")

setwd("~/Documents/Git-Projects/Git-Research-Projects/drug-response-prediction")
source("cnProfileVisualizationLibrary.R")
source("helperFunctions.R")

sample <- "hT30"

cd_doc()    
facets_segment_data <- retrieveFacetsSegments(sample, dir = "CSHL/Project_TUV_12995_B01_SOM_Targeted.2018-03-02/")
facets_segment_data <- segmentsToBedFormat(facets_segment_data)  

facets_snp_data <- retrieveFacetsSnps(sample, dir = "CSHL/Project_TUV_12995_B01_SOM_Targeted.2018-03-02/")
facets_snp_data <- snpsToBedFormat(facets_snp_data)

cd_local()
facets_segment_clusters <- retrieveFacetsSegmentClusters(sample, dir = "segClusteringResults/prev_run1/")
facets_segment_clusters <- segmentClustersToBedFormat(facets_segment_clusters)

visualizeCNProfile(facets_segment_data = facets_segment_data, facets_snp_data = facets_snp_data)

visualizeCNProfile(facets_segment_data = facets_segment_data, facets_snp_data = facets_snp_data, categories = "chr19")

visualizeCNProfile(facets_segment_data = facets_segment_clusters, facets_snp_data = facets_snp_data, categories = "chr19")




