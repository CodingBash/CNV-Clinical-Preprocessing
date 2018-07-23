
source("https://bioconductor.org/biocLite.R")
biocLite("gtrellis")
biocLite("ComplexHeatmap")

setwd("~/Documents/Git-Projects/Git-Research-Projects/drug-response-prediction")
source("cnProfileVisualizationLibrary.R")
source("helperFunctions.R")

sample <- "hN30"

setwd("~/Documents")
facets_segment_data_raw <- retrieveFacetsSegments(sample, dir = "CSHL/Project_TUV_12995_B01_SOM_Targeted.2018-03-02/")
facets_segment_data <- segmentsToBedFormat(facets_segment_data_raw)  
facets_segment_clusters_data <- segmentsToBedFormat(facets_segment_data_raw, median.clust = TRUE)

facets_snp_data_raw <- retrieveFacetsSnps(sample, dir = "CSHL/Project_TUV_12995_B01_SOM_Targeted.2018-03-02/")
facets_snp_data <- snpsToBedFormat(facets_snp_data_raw)

setwd("~/Documents/Git-Projects/Git-Research-Projects/drug-response-prediction")
facets_segment_clusters <- retrieveFacetsSegmentClusters(sample, dir = "segClusteringResults/prev_run1/")
facets_segment_clusters <- segmentClustersToBedFormat(facets_segment_clusters)

sample <- "hN30"
setwd("~/Documents")
facets_segment_data_raw <- retrieveFacetsSegments(sample, dir = "CSHL/Project_TUV_12995_B01_SOM_Targeted.2018-03-02/")
facets_segment_data <- segmentsToBedFormat(facets_segment_data_raw)
facets_snp_data_raw <- retrieveFacetsSnps(sample, dir = "CSHL/Project_TUV_12995_B01_SOM_Targeted.2018-03-02/")
facets_snp_data <- snpsToBedFormat(facets_snp_data_raw)
visualizeCNProfile(facets_segment_data = facets_segment_data, facets_snp_data = facets_snp_data)

sample <- "hN30"
setwd("~/Documents")
facets_snp_data_raw <- retrieveFacetsSnps(sample, dir = "CSHL/Project_TUV_12995_B01_SOM_Targeted.2018-03-02/")
facets_snp_data <- snpsToBedFormat(facets_snp_data_raw)
facets_snp_data_raw <- retrieveFacetsSnps(sample, dir = "CSHL/Project_TUV_12995_B01_SOM_Targeted.2018-03-02/")
facets_snp_data <- snpsToBedFormat(facets_snp_data_raw)
visualizeCNProfile(facets_segment_data = facets_segment_data, facets_snp_data = facets_snp_data, categories = "chr19")

sample <- "hN30"
setwd("~/Documents")
facets_segment_data_raw <- retrieveFacetsSegments(sample, dir = "CSHL/Project_TUV_12995_B01_SOM_Targeted.2018-03-02/")
facets_segment_clusters_data <- segmentsToBedFormat(facets_segment_data_raw, median.clust = TRUE)
facets_snp_data_raw <- retrieveFacetsSnps(sample, dir = "CSHL/Project_TUV_12995_B01_SOM_Targeted.2018-03-02/")
facets_snp_data <- snpsToBedFormat(facets_snp_data_raw)
visualizeCNProfile(facets_segment_data = facets_segment_clusters_data, facets_snp_data = facets_snp_data, categories = "chr19")

sample <- "hT30"
setwd("~/Documents")
facets_segment_data_raw <- retrieveFacetsSegments(sample, dir = "CSHL/Project_TUV_12995_B01_SOM_Targeted.2018-03-02/")
facets_segment_data <- segmentsToBedFormat(facets_segment_data_raw)
facets_snp_data_raw <- retrieveFacetsSnps(sample, dir = "CSHL/Project_TUV_12995_B01_SOM_Targeted.2018-03-02/")
facets_snp_data <- snpsToBedFormat(facets_snp_data_raw)
visualizeCNProfile(facets_segment_data = facets_segment_data, facets_snp_data = facets_snp_data)



sample <- "hT30"
setwd("~/Documents")
facets_segment_data_raw <- retrieveFacetsSegments(sample, dir = "CSHL/Project_TUV_12995_B01_SOM_Targeted.2018-03-02/")
facets_segment_clusters_data <- segmentsToBedFormat(facets_segment_data_raw, median.clust = TRUE)
facets_snp_data_raw <- retrieveFacetsSnps(sample, dir = "CSHL/Project_TUV_12995_B01_SOM_Targeted.2018-03-02/")
facets_snp_data <- snpsToBedFormat(facets_snp_data_raw)
visualizeCNProfile(facets_segment_data = facets_segment_clusters_data, facets_snp_data = facets_snp_data)



sample <- "hT30"
setwd("~/Documents/Git-Projects/Git-Research-Projects/drug-response-prediction")
facets_segment_clusters_raw <- retrieveFacetsSegmentClusters(sample, dir = "segClusteringResults/prev_run1/")
facets_segment_clusters <- segmentClustersToBedFormat(facets_segment_clusters_raw)
setwd("~/Documents")
facets_snp_data_raw <- retrieveFacetsSnps(sample, dir = "CSHL/Project_TUV_12995_B01_SOM_Targeted.2018-03-02/")
facets_snp_data <- snpsToBedFormat(facets_snp_data_raw)
visualizeCNProfile(facets_segment_data = facets_segment_clusters, facets_snp_data = facets_snp_data)

sample <- "hT58"
setwd("~/Documents")
facets_segment_data_raw <- retrieveFacetsSegments(sample, dir = "CSHL/Project_TUV_12995_B01_SOM_Targeted.2018-03-02/")
facets_segment_data <- segmentsToBedFormat(facets_segment_data_raw)
facets_snp_data_raw <- retrieveFacetsSnps(sample, dir = "CSHL/Project_TUV_12995_B01_SOM_Targeted.2018-03-02/")
facets_snp_data <- snpsToBedFormat(facets_snp_data_raw)
visualizeCNProfile(facets_segment_data = facets_segment_data)


