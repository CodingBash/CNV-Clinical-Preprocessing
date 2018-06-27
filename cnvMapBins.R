
cd_local <- function() {
  current_path <- getActiveDocumentContext()$path 
  setwd(dirname(current_path ))
}

cd_doc <- function() {
  setwd("C:/Users/bbece/Documents")
}

cd_local()
samples <- read.table("sampleList.csv", header=T, sep = "\t", stringsAsFactors = F)
#for(sample in samples$Organoids){
sample <- "hN31"

# Import data
cd_doc()

facets_data <- as.data.frame(read.table(paste("CSHL/Project_TUV_12995_B01_SOM_Targeted.2018-03-02/Sample_", sample, "/analysis/structural_variants/", sample, "--NA12878.cnv.facets.v0.5.2.txt", sep = ""), header = TRUE, sep="\t", stringsAsFactors=FALSE, quote=""))
facets_data <- facets_data[,c(1, 10, 11, 5)]
names(facets_data) <- c("chr", "start", "end", "cnlr")
facets_data[facets_data$chr == "X", ]$chr <- "23" #TODO: May not work if multiple segments in X chromosome
facets_data[facets_data$chr == "Y", ]$chr <- "24"
facets_data$chr <- as.numeric(facets_data$chr)

facets_bins_data <- as.data.frame(read.table(paste("CSHL/Project_TUV_12995_B01_SOM_Targeted.2018-03-02/Sample_", sample, "/analysis/structural_variants/", sample, "--NA12878.procSample-jseg.cnv.facets.v0.5.2.txt", sep = ""), header = TRUE, sep="\t", stringsAsFactors=FALSE, quote=""))
facets_bins_data$binNumber <- 1:nrow(facets_bins_data)
facets_bins_data <- facets_bins_data[,c(1, 2, 17)]
names(facets_bins_data) <- c("chr", "maploc", "binNumber")
facets_bins_data$chr <- as.numeric(facets_bins_data$chr)

for(segment in 1:nrow(facets_data)){
  start_binNumber <- facets_bins_data[facets_bins_data$chr == facets_data[segment, ]$chr & facets_bins_data$maploc == facets_data[segment,]$start, ]$binNumber + 1 # Finds bin number with start == maploc, then gets the true bin after
  end_binNumber <- facets_bins_data[facets_bins_data$chr == facets_data[segment, ]$chr & facets_bins_data$maploc == facets_data[segment,]$end, ]$binNumber + 1 # Finds bin number with end == maploc, then gets the true bin after 
  facets_data[segment,]$start <- start_binNumber
  facets_data[segment,]$end <- end_binNumber
}
