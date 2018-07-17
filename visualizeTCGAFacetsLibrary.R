library(facets)

plotFacetsSample <- function(tumorId, normalId, fitPrefix = "facetsG5Fit_", ooPrefix = "facetsG5OO_", silent = TRUE){
  try({
    fit <- readRDS(paste(fitPrefix, tumorId, "_", normalId, ".rds", sep = ""))
    oo <- readRDS(paste(ooPrefix, tumorId, "_", normalId, ".rds", sep = ""))
    plotSample(x=oo, emfit = fit)
  }, silent = silent)
}