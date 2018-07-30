#
# This script handles the conversion of bases to bins #mapBasesToBins() and bins to bases (using a binning scheme)
#
setwd("~/Git-Projects/Git-Research-Projects/drug-response-prediction")
source("helperFunctions.R")

reference <- "hN31"
sample_dir <- "./output/FACETS_Reference_hN31_7_28_18_2/"


#
# Convert bp to bins from a set of sample files
#
mapBasesToBins <- function(){
  #
  # Map samples
  #
  mapSample <- function(sample){
    setwd("~/Git-Projects/Git-Research-Projects/FACETS_write_files")
    # TODO: prevent duplicate code by using #retrieveBinningScheme()
    
    #
    # Preprocess segments data
    #
    facets_data <- retrieveFacetsSegments(sample, sample_subdir = "/", reference = "hN31", dir = sample_dir)
    facets_data <- facets_data[,c(1, 10, 11, 5)]
    names(facets_data) <- c("chr", "start", "end", "cnlr")

    if (length(facets_data[facets_data$chr == "X", ]$chr) > 0) facets_data[facets_data$chr == "X", ]$chr <- "23" #TODO: May not work if multiple segments in X chromosome
    if (length(facets_data[facets_data$chr == "Y", ]$chr) > 0)facets_data[facets_data$chr == "Y", ]$chr <- "24"
    facets_data$chr <- as.numeric(facets_data$chr)
    
    print(paste0("Preprocess segments data for sample=", sample))
    #
    # Preprocess SNPs data
    #
    facets_bins_data <- retrieveFacetsSnps(sample, sample_subdir = "/", reference = "hN31", dir = sample_dir)
    facets_bins_data$binNumber <- 1:nrow(facets_bins_data)
    facets_bins_data <- facets_bins_data[,c(1, 2, ncol(facets_bins_data))]
    names(facets_bins_data) <- c("chr", "maploc", "binNumber")
    facets_bins_data$chr <- as.numeric(facets_bins_data$chr)
    print(paste0("Preprocess SNPs data for sample=", sample))
    #
    # Assign bin number
    #
    for(segment in 1:nrow(facets_data)){
      print(segment)
      start_binNumber <- facets_bins_data[facets_bins_data$chr == facets_data[segment, ]$chr & facets_bins_data$maploc == facets_data[segment,]$start, ]$binNumber + 1 # Finds bin number with start == maploc, then gets the true bin after
      end_binNumber <- facets_bins_data[facets_bins_data$chr == facets_data[segment, ]$chr & facets_bins_data$maploc == facets_data[segment,]$end, ]$binNumber + 1 # Finds bin number with end == maploc, then gets the true bin after 
      print(end_binNumber)
      facets_data[segment,]$start <- start_binNumber
      facets_data[segment,]$end <- end_binNumber
    }
    
    print(paste("Mapped ", sample, sep = ""))
    return(facets_data)
    
  }
  
  outputMappedSample <- function(sample, mapped_sample){
    cd_local()
    file_name <- paste("mappedFacetsFiles/", sample, "--", reference, ".mapped.cnv.facets.v0.5.2.bed", sep = "")
    write.csv(mapped_sample, file = file_name, row.names = FALSE)
    print(paste("Wrote ", sample, " to ", file_name))
  }
  
  #
  # Read in samples to map
  #
  cd_local()
  loaded_samples <- load_samples(classes = c("T", "F", "M"), sampleList = "sampleList.csv")
    
  for(sample in loaded_samples){
    mapped_sample <- mapSample(sample)
    print(head(mapped_sample))
    outputMappedSample(sample, mapped_sample)
  }
}


retrieveBinningScheme = function(sample = "hN31") {
  cd_doc()
  facets_bins_data <- as.data.frame(read.table(paste("CSHL/Project_TUV_12995_B01_SOM_Targeted.2018-03-02/Sample_", sample, "/analysis/structural_variants/", sample, "--NA12878.procSample-jseg.cnv.facets.v0.5.2.txt", sep = ""), header = TRUE, sep="\t", stringsAsFactors=FALSE, quote=""))
  facets_bins_data$binNumber <- 1:nrow(facets_bins_data)
  facets_bins_data <- facets_bins_data[,c(1, 2, 17)]
  names(facets_bins_data) <- c("chr", "maploc", "binNumber")
  facets_bins_data$chr <- as.numeric(facets_bins_data$chr)
  return(facets_bins_data)
}

#
# For a dataframe from a single sample, map bins to bases
# - Knowing sample is potentially important since binning scheme may different
# - ^ This may be a challenge down the line since sample may be unknown in a concatenated datafram
# TODO: Do this once it is time for visualization
#
#
# @param bedDataFrame - dataframe of intervals in bed format (chrom (atomic vector), start (bins), end (bins))
# @param binningScheme - dataframe of binning scheme (bin number, chr, bp)
#
mapBinsToBases <- function(bedDataframe, binningScheme = retrieveBinningScheme()) {
  for(interval in 1:nrow(bedDataframe)){
    start_baseNumber <- binningScheme[binningScheme$binNumber == as.integer(bedDataframe[interval, ]["start"]) - 1, ]$maploc
    end_baseNumber <- binningScheme[binningScheme$binNumber == as.integer(bedDataframe[interval, ]["end"]) - 1, ]$maploc
    bedDataframe[interval, ]["start"] <- start_baseNumber
    bedDataframe[interval, ]["end"] <- end_baseNumber
  }
  return(bedDataframe)
}

# mapBasesToBins()
