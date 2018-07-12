args <- commandArgs(trailingOnly = TRUE)

#
# Default script arguments
#
event <- "A" 
outputCsv <- "coreTable.csv"
outputObj <- "newCOREobj.rds"
if (length(args) == 1){
	event <- args[1]
} else if (length(args) == 2){
	event <- args[1]
	outputCsv <- args[2]	
} else if (length(args) == 3){
	event <- args[1]
	outputCsv <- args[2]
	outputObj <- args[3]
}

# TODO: These may not work.
source("coreGenerationLibrary.R")
source("helperFunctions.R")

cd_home()
samples <- load_samples(classes = c("N"), sampleList = "./resources/sampleList.csv")
chromosomeSizes <- readRDS("./resources/chromosomeSizes.rds")
inputCORESegments <- generateInputCORESegments(event, samples, chromosomeSizes, dir = "./resources/Project_TUV_12995_B01_SOM_Targeted.2018-03-02/", extension = "cnv.facets.v0.5.2.txt", rescaleInput = TRUE, ampCall = 0.2, delCall = -0.235)
inputCOREBoundaries <- generateInputCOREBoundaries(chromosomeSizes)
print("Prepared all inputs - now running CORE")

outputCOREobj <- runCORE(inputCORESegments, inputCOREBoundaries, distrib="Grid", maxmark=500, nshuffle=500, seedme, njobs=4)
print("CORE run complete")

saveRDS(outputCOREobj, outputObj)
print(paste("Saved CORE obj to", outputObj, sep = ""))

retrieveCORETable(outputCOREobj, chromosomeSizes, rescaleOutput = TRUE)
write.csv(coreTable, outputCsv)
print(paste("Saved coreTable as", outputCsv, sep = ""))

print("script complete")
