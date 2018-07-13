#
# Quick interactive script to sub-select cores based on p-value
# TODO: Make into callable script with arguments
#

#
# Load source libraries
#
setwd("~/Git-Projects/Git-Research-Projects/drug-response-prediction")
source("helperFunctions.R")

#
# Set the core p_threshold
#
p_threshold <- 0.002

#
# Get core information
#
cd_local()
Aobj <- readRDS("./hN_results/coresResults/AnewCOREobj.rds") # Retrieve CORE object
coreTable <- read.table("AcoreTableBP.csv", header = TRUE, sep = ",") # Retrieve CORE table CSV (since the table's scale is chromosome-based instead of absolute)

#
# Subset core information using p_threshold then write to file
#
coreTable <- coreTable[which(Aobj$p<p_threshold),] # Filter cores based on p-value threshold
coreTable <- coreTable[,c(2,3,4)] # Now convert to BED format
write.table(coreTable, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE, file = "AcoresBP.bed") # Write to bed file