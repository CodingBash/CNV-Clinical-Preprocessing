
setwd("~/Documents/Git-Projects/Git-Research-Projects/drug-response-prediction")
library(facets)
source("facetsAnalysisLibrary.R")

plotFacetsSample <- function(tumorId, normalId, silent = TRUE){
  try({
    setwd(paste0("~/Documents/Git-Projects/Git-Research-Projects/FACETS_nonmatching_test/"))
    fit <- getFacetsFit(tumorId, normalId, dir = "output/prev_run_4/")
    oo <- getFacetsOO(tumorId, normalId, dir = "output/prev_run_4/")
    print(paste0("FACETS Profile for tumorId=", tumorId, " normalId=", normalId))
    plotSample(x=oo, emfit = fit)
  }, silent = silent)
}
      

plotFacetsSample(1, 1)

plotFacetsSample(1, 2)

plotFacetsSample(1, 5)

plotFacetsSample(2, 2)

plotFacetsSample(2, 1)

plotFacetsSample(2, 3)

plotFacetsSample(2, 4)

plotFacetsSample(5, 5)

plotFacetsSample(5, 6)

plotFacetsSample(5, 7)

plotFacetsSample(6, 6)

plotFacetsSample(6, 7)

plotFacetsSample(6, 5)

plotFacetsSample(7, 7)

plotFacetsSample(7, 1)

plotFacetsSample(7, 2)


