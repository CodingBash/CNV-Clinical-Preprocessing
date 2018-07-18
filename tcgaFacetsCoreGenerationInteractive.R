setwd("~/Git-Projects/Git-Research-Projects/drug-response-prediction")
source("helperFunctions.R")
source("coreGenerationLibrary.R")
source("facetsAnalysisLibrary.R")
library(BSgenome.Hsapiens.UCSC.hg19)

#
# Calculating CORES of non matching FACET pairs
#
for(tumorId in seq(1, 7)) {
  cd_facets()
  fitFiles <- getAllFitFilesForTumor(tumorId, dir = "output/prev_run_1/", bedFormat = FALSE)
  chromosomeSizes <- generateChromosomeSizes(BSgenome.Hsapiens.UCSC.hg19)
  events = c("A")
  inputCORESegments <- subsetAllSegmentsByEvent(fitFiles, events, chromosomeSizes, rescaleInput = TRUE, ampCall = 0.2, delCall = -0.235)
  
  inputCOREBoundaries <- generateInputCOREBoundaries(chromosomeSizes)
  
  #
  # Run CORE
  #
  outputCOREobj <- runCORE(inputCORESegments, inputCOREBoundaries, distrib="Rparallel", maxmark=10, nshuffle=20, seedme=123, njobs=4)
  
  #
  # Get CORE table
  #
  COREtable <- retrieveCORETable(outputCOREobj, chromosomeSizes, rescaleOutput = TRUE)
  
}
