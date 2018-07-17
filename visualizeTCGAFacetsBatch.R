setwd("~/Git-Projects/Git-Research-Projects/drug-response-prediction")
source("visualizeTCGAFacetsLibrary.R")
source("helperFunctions.R")

for(tumorId in seq(1, 7)){
  for(normalId in seq(1, 7)){
    cd_local()
    pdf(paste0("tcgaFacetsVisualization/tcgaFacets_t", tumorId, "_n", normalId, ".pdf"),width=16,height=9) 
    cd_facets("output/")
    plotFacetsSample(tumorId, normalId)
    dev.off()
  }
}

