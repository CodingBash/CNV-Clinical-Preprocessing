
source("https://bioconductor.org/biocLite.R")
biocLite("gtrellis")
biocLite("ComplexHeatmap")

source("cnProfileVisualizationLibrary.R")
source("helperFunctions.R")

sample <- "hT30"

cd_doc()    
facets_segment_data <- retrieveFacetsSegments(sample, dir = "CSHL/Project_TUV_12995_B01_SOM_Targeted.2018-03-02/")
facets_snp_data <- retrieveFacetsSnps(sample, dir = "CSHL/Project_TUV_12995_B01_SOM_Targeted.2018-03-02/")
facets_segment_data <- segmentsToBedFormat(facets_segment_data)  
facets_snp_data <- snpsToBedFormat(facets_snp_data)

visualizeCNProfile(facets_segment_data = facets_segment_data, facets_snp_data = facets_snp_data)

visualizeCNProfile(facets_segment_data = facets_segment_data, facets_snp_data = facets_snp_data, categories = "chr19")


