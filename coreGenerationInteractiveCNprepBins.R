#
# This script generates CORES from input files in BP units (instead of SNP/bin units).
# This script is interactive and is NOT meant to be ran on a HPC job from a shell script.
#

#
# Only run if not installed
#
#source("https://bioconductor.org/biocLite.R")
#biocLite("BSgenome.Hsapiens.UCSC.hg19")
#install.packages("CORE")

library(BSgenome.Hsapiens.UCSC.hg19)

setwd("~/Git-Projects/Git-Research-Projects/drug-response-prediction")
source("coreGenerationLibrary.R")
source("helperFunctions.R")

events <- c("D") # "A" - amplification, "D" - deletion


#
# Get CORE input
#
cd_local()
samples <- load_samples(classes = c("T"), sampleList = "sampleList.csv")

setwd("~/Git-Projects/Git-Research-Projects/cnprep_clustering")

inputCORESegments <- selectCnprepSegmentsWithEvent(events, samples, "output/prev_run_7_27_2018_8", 0.001, probes = TRUE, silent = FALSE)

#
# Run CORE
#
outputCOREobj <- runCORE(inputCORESegments, distrib="Rparallel", maxmark=10, nshuffle=20, seedme=123, njobs=4)

#
# Get CORE table
#
COREtable <- retrieveCORETable(outputCOREobj, chromosomeSizes, rescaleOutput = TRUE)
