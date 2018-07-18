#
# This script generates CORES for the TCGA FACETS segments across all nonmatching tumor-normal profiles
# for each tumor ID. Therefore, we can compare the COREs to the matching tumor-normal to see if they are similar or not.
#
args <- commandArgs(trailingOnly = TRUE)

#
# Set script arguments
#
event <- "A" 
outputCsv <- "coreTable" #_t*.csv appended
outputObj <- "newCOREobj" #_t*.rds appended
tumorId <- NA
if (length(args) == 1){
  event <- args[1]
} else if (length(args) == 2){
  event <- args[1]
  outputCsv <- args[2]	
} else if (length(args) == 3){
  event <- args[1]
  outputCsv <- args[2]
  outputObj <- args[3]
}  else if (length(args) == 4){
  event <- args[1]
  outputCsv <- args[2]
  outputObj <- args[3]
  tumorId <- as.numeric(args[4])
}

source("helperFunctions.R")
source("coreGenerationLibrary.R")
source("facetsAnalysisLibrary.R")

#
# Calculating CORES of non matching FACET pairs
#
chromosomeSizes <- readRDS("./resources/chromosomeSizes.rds")
fitFiles <- getAllFitFilesForTumor(tumorId, dir = "./", bedFormat = FALSE)
events = c("A")
inputCORESegments <- subsetAllSegmentsByEvent(fitFiles, events, chromosomeSizes, rescaleInput = TRUE, ampCall = 0.2, delCall = -0.235)
inputCOREBoundaries <- generateInputCOREBoundaries(chromosomeSizes)
print("Prepared all inputs - now running CORE")

#
# Run CORE
#
outputCOREobj <- runCORE(inputCORESegments, inputCOREBoundaries, distrib="Grid", maxmark=500, nshuffle=500, seedme=123, njobs=4)
print("CORE run complete")

#
# Save CORE object
# WARNING: The coreTable in the CORE object may not be in chromosomal location units. See CORE table output
#
saveRDS(outputCOREobj, paste0(outputObj, "_t", tumorId, ".rds"))
print(paste("Saved CORE obj to", paste0(outputObj, "_t", tumorId, ".rds"), sep = " "))

#
# Get CORE table
#
COREtable <- retrieveCORETable(outputCOREobj, chromosomeSizes, rescaleOutput = TRUE)
write.csv(CORETable, paste0(outputCsv, "_t", tumorId, ".csv"))
print(paste("Saved coreTable as", paste0(outputCsv, "_t", tumorId, ".csv"), sep = " "))

print("script complete")
