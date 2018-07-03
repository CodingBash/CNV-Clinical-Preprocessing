Aobj <- readRDS("C:/Users/bbece/Documents/Git-Projects/Git-Research-Projects/drug-response-prediction/hN_results/coresResults/AnewCOREobj.rds")

coreTable <- read.table("AcoreTableBP.csv", header = TRUE, sep = ",")
coreTable <- coreTable[which(Aobj$p<0.002),]
coreTable <- coreTable[,c(2,3,4)] # Now convert to BED format
write.table(coreTable, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE, file = "AcoresBP.bed")