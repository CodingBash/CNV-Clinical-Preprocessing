#
# Select Tcga FACET cores based on p value
#
setwd("~/Git-Projects/Git-Research-Projects/drug-response-prediction")
source("helperFunctions.R")


selectCores <- function(event, tumorId){
  try({
    cd_tcga("output")
    Aobj <- readRDS(paste0(event, "tcgaCOREobjBP_t", tumorId, ".rds"))
    coreTable <- read.table(paste0(event, "tcgaCoreTableBP_t", tumorId, ".csv"), header = TRUE, sep = ",")
    coreTable <- coreTable[which(Aobj$p<0.002),]
    coreTable <- coreTable[,c(2,3,4)] # Now convert to BED format
    
    cd_facets("cores")
    write.table(coreTable, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE, file = paste0(event, "coresBP_t", tumorId, ".bed"))    
  }, silent = TRUE)
}
events <- c("A", "D")
tumorIds <- seq(1,7)
for(event in events){
  for(tumorId in tumorIds){
    selectCores(event, tumorId)
  }
}


