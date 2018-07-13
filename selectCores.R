#
# Quick interactive script to sub-select cores based on p-value
# TODO: Make into callable script with arguments
#
source("helperFunctions.R")

p_threshold <- 0.002

cd_local()
Aobj <- readRDS("./hN_results/coresResults/AnewCOREobj.rds") # Retrieve CORE object
coreTable <- read.table("AcoreTableBP.csv", header = TRUE, sep = ",") # Retrieve CORE table CSV (since the table's scale is chromosome-based instead of absolute)
coreTable <- coreTable[which(Aobj$p<p_threshold),] # Filter cores based on p-value threshold
coreTable <- coreTable[,c(2,3,4)] # Now convert to BED format
write.table(coreTable, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE, file = "AcoresBP.bed") # Write to bed file