#source("https://bioconductor.org/biocLite.R")
#biocLite("GenomicRanges")

setwd("~/Git-Projects/Git-Research-Projects/drug-response-prediction")
source("helperFunctions.R")
library(GenomicRanges)

#
# Load sample to retrieve feature set for
#
samples <- load_samples(classes = c("T"))
for(sample in samples){
  
  #
  # Retrieve feature labels (COREs) analgous to all samples
  cd_local("hT_results/finalCores/prev_run_7_3_2018")
  Acores <- read.table(paste0("AfinalCoresBP_subtractFraction.bed"), header = FALSE, sep = "\t", stringsAsFactors = FALSE) 
  Dcores <- read.table(paste0("DfinalCoresBP_subtractFraction.bed"), header = FALSE, sep = "\t", stringsAsFactors = FALSE) 
  
  #
  # Retrieve necessary data for feature value calculation (bins with cnlr and COREs)
  #
  cd_doc()
  facets_snp_data <- retrieveFacetsSnps(sample, dir = "CSHL/Project_TUV_12995_B01_SOM_Targeted.2018-03-02/")
  snp_bed <- snpsToBedFormat(facets_snp_data)
  
  
  #
  # Preprocess COREs
  #
  seqnames_A <- table(Acores[[1]])
  gr_A <- GRanges(
    seqnames = Rle(names(seqnames_A), as.vector(seqnames_A)),
    ranges = IRanges(Acores[[2]], Acores[[3]], names = row.names(Acores)),
    event=rep("A", nrow(Acores))
  )
  seqnames_D <- table(Dcores[[1]])
  gr_D <- GRanges(
    seqnames = Rle(names(seqnames_D), as.vector(seqnames_D)),
    ranges = IRanges(Dcores[[2]], Dcores[[3]], names = row.names(Dcores)),
    event=rep("D", nrow(Dcores))
  )
  # TODO: Do operations on GRanges object to simplify, if need be (i.e. reduction)
  gr <- c(gr_A, gr_D)
  coreDf <- as.data.frame(gr, row.names = seq_along(gr$event))
  coreDf$cnlr <- NA
  
  #
  # Determine value of each CORE feature (cnlr)
  #
  coreDf[, 7] <- sapply(seq(1, nrow(coreDf)), function(core.index){
    # Get snps with same chromosome as core, and in between the core region
    snps <- snp_bed[snp_bed[[1]] == coreDf[core.index, 1] & snp_bed[[2]] >= coreDf[core.index, 2] & snp_bed[[3]] <= coreDf[core.index, 3], ]
    # Calculate and assign median
    return(median(snps$value))
  } )
  
  # Convert back to Granges object - this is our final feature object
  featureSet <- makeGRangesFromDataFrame(coreDf, keep.extra.columns = TRUE)
  print(featureSet)
}
