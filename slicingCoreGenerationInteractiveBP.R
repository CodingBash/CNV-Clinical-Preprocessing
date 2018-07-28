#
# This script generates CORES from input files in BP units (instead of SNP/bin units).
# This script is interactive and is NOT meant to be ran on a HPC job from a shell script.
#

#
# Only run if not installed
#
source("https://bioconductor.org/biocLite.R")
biocLite("BSgenome.Hsapiens.UCSC.hg19")
install.packages("CORE")

library(BSgenome.Hsapiens.UCSC.hg19)

setwd("~/Git-Projects/Git-Research-Projects/drug-response-prediction")
source("coreGenerationLibrary.R")
source("helperFunctions.R")

events <- c("A", "D") # "A" - amplification, "D" - deletion

#
# Get CORE input
#
cd_local()
samples <- load_samples(classes = c("T", "M", "F"), sampleList = "sampleList.csv")
chromosomeSizes <- generateChromosomeSizes(genome = BSgenome.Hsapiens.UCSC.hg19)

cd_local()
inputCORESegments <- loadSlicingRegions(dir = "slicingOutput/prev_run_7_28_18_4/", samples = samples, events = events, silent = TRUE)
inputCOREBoundaries <- generateInputCOREBoundaries(chromosomeSizes)

#
# Run CORE
#
outputCOREobj <- runCORE(inputCORESegments, inputCOREBoundaries, distrib="Rparallel", maxmark=10, nshuffle=20, seedme=123, njobs=4)

#
# Get CORE table
#
COREtable <- retrieveCORETable(outputCOREobj, chromosomeSizes, rescaleOutput = TRUE)
