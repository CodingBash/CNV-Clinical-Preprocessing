#
# Get FACETS xx file
#
getFacetsXX <- function(tumorId, normalId, xxPrefix = "facetsG5XX_", dir = "output/"){
  xx <- readRDS(paste0(dir, xxPrefix, tumorId, "_", normalId, ".rds"))
  return(xx)
}

#
# Get FACETS oo file
#
getFacetsOO <- function(tumorId, normalId, ooPrefix = "facetsG5OO_", dir = "output/"){
  oo <- readRDS(paste0(dir, ooPrefix, tumorId, "_", normalId, ".rds"))
  return(oo)
}

#
# Get FACETS fit file
#
getFacetsFit <- function(tumorId, normalId, fitPrefix = "facetsG5Fit_", dir = "output/"){
  fit <- readRDS(paste0(dir, fitPrefix, tumorId, "_", normalId, ".rds"))
  return(fit)
}

# TODO: Redundant function to get CNCF attribute from fit object. 
retrieveFacetsSegmentsFromObject <- function(tumorId, normalId, fitPrefix = "facetsG5Fit_", dir = "output/"){
  fit <- readRDS(paste(dir, fitPrefix, tumorId, "_", normalId, ".rds", sep = ""))
  return(fit$cncf)
}

#
# Get all SNP profiles as list for a specific tumor ID and all normal IDs (matching or nonmatching)
#
getAllXXFilesForTumor <- function(tumorId, xxPrefix = "facetsG5XX_", dir = "output/prev_run_1/", bedFormat = TRUE, silent = FALSE){
  xxList <- list()
  for(normal.id in seq(1,7)){
    try({
      xx <- getFacetsXX(tumorId, normal.id, xxPrefix = xxPrefix, dir = dir)
      if(bedFormat == TRUE){
        snps <- snpsToBedFormat(xx$jointseg)
        xxList[[normal.id]] <- snps
      } else {
        xxList[[normal.id]] <- xx$jointseg
      }
    }, silent = silent)
  }
  return(xxList)
}

#
# Get all segment profiles as list for a specific tumor ID and all normal IDs (matching or nonmatching)
#
getAllFitFilesForTumor <- function(tumorId, fitPrefix = "facetsG5Fit_", dir = "output/prev_run_1/", bedFormat = TRUE, silent = FALSE){
  fitList <- list()
  for(normal.id in seq(1,7)){
    try({
    fit <- getFacetsFit(tumorId, normal.id, fitPrefix = fitPrefix, dir = dir)
      if(bedFormat == TRUE){
        segments <- segmentsToBedFormat(fit$cncf)
        fitList[[normal.id]] <- segments
      } else {
        fitList[[normal.id]] <- fit$cncf
      }
    }, silent = silent)
  }
  return(fitList)
}
