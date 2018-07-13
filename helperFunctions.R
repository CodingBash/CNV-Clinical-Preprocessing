#
# This script contains a set of helper functions for the other CNA scripts
#

#
# Sets the directory to the project workspace
# TODO: directory is hardcoded
#
cd_local <- function() {
  setwd("C:/Users/bbece/Documents/Git-Projects/Git-Research-Projects/drug-response-prediction")
}

#
# Sets the directory to the project workspace
# TODO: directory is hardcoded
#
cd_doc <- function() {
  setwd("C:/Users/bbece/Documents")
}

#
# Sets the directory to HPC home directory
# TODO: migrate off of this function and use cd_core() instead
#
cd_home <- function() {
  setwd("~/code/hN_core_artifacts")
}

cd_core <- function() {
  setwd("~/code/hN_core_artifacts")  
}

cd_cnprep <- function() {
  setwd("~/code/cnprep_clustering")
}

# 
# Load all sample names given a sample class
#
load_samples <- function(classes = c("N"), sampleList = "sampleList.csv") {
  samples <- read.table(sampleList, header=T, sep = "\t", stringsAsFactors = F)
  loaded_samples <- c(NA)
  loaded_samples.index <- 1
  for(sample in samples$Organoids){
    for(class in classes){
      if(substring(sample, 2,2) == class){
        loaded_samples[loaded_samples.index] <- sample
        loaded_samples.index <- loaded_samples.index + 1
        next
      }
    }
  }
  return(loaded_samples)
}

#
# Retrieve a sample's FACETS segmented data from specified directory
#
retrieveFacetsSegments <- function(sample, dir = "CSHL/Project_TUV_12995_B01_SOM_Targeted.2018-03-02/"){
  facets_segment_data <- as.data.frame(read.table(paste(dir, "Sample_", sample, "/analysis/structural_variants/", sample, "--NA12878.cnv.facets.v0.5.2.txt", sep = ""), header = TRUE, sep="\t", stringsAsFactors=FALSE, quote=""))
  return(facets_segment_data)
}

#
# Retrieve a sample's FACETS segmented cluster (cnlr.median.clust) data from specified directory
#
retrieveFacetsSegmentClusters <- function(sample, dir = "segClusteringResults/"){
  facets_segment_clusters <- as.data.frame(read.table(paste(dir, sample, "_segtable.tsv", sep = ""), header = TRUE, sep="\t", stringsAsFactors=FALSE, quote=""))
  return(facets_segment_clusters)
}

#
# Retrieve a sample's FACETS SNP data from specified directory
#
retrieveFacetsSnps <- function(sample, dir = "CSHL/Project_TUV_12995_B01_SOM_Targeted.2018-03-02/"){
  facets_snp_data <- as.data.frame(read.table(paste(dir, "Sample_", sample, "/analysis/structural_variants/", sample, "--NA12878.procSample-jseg.cnv.facets.v0.5.2.txt", sep = ""), header = TRUE, sep="\t", stringsAsFactors=FALSE, quote=""))
  return(facets_snp_data)
}

#
# Simplifies FACETS segment original format into BED format with
# the columns: "chrom", "chrom.start", "chrom.end", "median cnlr"
#
segmentsToBedFormat <- function(facets_segment_data){
  facets_segment_data <- facets_segment_data[,c(1, 10, 11, 5)]
  facets_segment_data$X.chrom. <- paste("chr", facets_segment_data$X.chrom., sep="")
  names(facets_segment_data) <- c("chrom", "start", "end", "value")
  return(facets_segment_data)
}

#
# Simplifies FACETS segment cluster original format into BED format with
# the columns: "chrom", "chrom.start", "chrom.end", "median cnlr"
#
segmentClustersToBedFormat <- function(facets_segment_clusters){
  facets_segment_clusters <- facets_segment_clusters[,c(6,7,8,20)]
  facets_segment_clusters <- facets_segment_clusters[facets_segment_clusters$chrom != "X",]
  facets_segment_clusters$chrom <- paste("chr", facets_segment_clusters$chrom, sep = "")
  names(facets_segment_clusters) <- c("chrom", "start", "end", "value")
  return(facets_segment_clusters)
}

#
# Simplifies FACETS snp original format into BED format. Also "bin-ifies" the data
# start = previous maploc, and end = current maploc
# Columns: "chrom", "new bin.start", "new bin.end", "cnlr"
#
snpsToBedFormat <- function(facets_snp_data){
  facets_snp_data$start <- seq(1, length.out=nrow(facets_snp_data), by=1)
  facets_snp_data$end <- seq(2, length.out=nrow(facets_snp_data), by=1)
  facets_snp_data <- facets_snp_data[,c(1, 2, 2, 11)]
  facets_snp_data$X.chrom. <- paste("chr", facets_snp_data$X.chrom., sep="")
  names(facets_snp_data) <- c("chrom", "start", "end", "value")
  return(facets_snp_data)
}
