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

events <- c("D") # "A" - amplification, "D" - deletion

sample_dir <- "output/FACETS_Reference_hN31_7_28_18_2/"

#
# Get CORE input
#
cd_local()
samples <- load_samples(classes = c("T"), sampleList = "sampleList.csv")
chromosomeSizes <- generateChromosomeSizes(genome = BSgenome.Hsapiens.UCSC.hg19)
setwd("~/Git-Projects/Git-Research-Projects/FACETS_write_files/")
#inputCORESegments <- selectSegmentsWithEvents(events = events, samples = samples, chromosomeSizes = chromosomeSizes, 
#                                             dir = "CSHL/Project_TUV_12995_B01_SOM_Targeted.2018-03-02/", extension = "cnv.facets.v0.5.2.txt", inSampleFolder = TRUE, 
#                                             rescaleInput = TRUE, ampCall = 0.2, delCall = -0.235)
inputCORESegments <- selectSegmentsWithEvents(events = events, samples = samples, chromosomeSizes = chromosomeSizes,
                                              dir = sample_dir, sample_subdir="/", reference=reference, extension = "cnv.facets.v0.5.2.txt", inSampleFolder = TRUE,
                                              rescaleInput = TRUE, ampCall = 0.2, delCall = -0.235)

inputCOREBoundaries <- generateInputCOREBoundaries(chromosomeSizes)

#
# Run CORE
#
outputCOREobj <- runCORE(inputCORESegments, inputCOREBoundaries, distrib="Rparallel", maxmark=10, nshuffle=20, seedme=123, njobs=4)

#
# Get CORE table
#
COREtable <- retrieveCORETable(outputCOREobj, chromosomeSizes, rescaleOutput = TRUE)
