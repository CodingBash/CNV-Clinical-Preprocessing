#
# Create fuctions to change directory
#
cd_local <- function() {
  current_path <- getActiveDocumentContext()$path 
  setwd(dirname(current_path ))
}

cd_doc <- function() {
  setwd("C:/Users/bbece/Documents")
}

#
# Convert bp to bins from a set of sample files
#
mapBasesToBins <- function(){
  #
  # Read in samples to map
  #
  cd_local()
  samples <- read.table("sampleList.csv", header=T, sep = "\t", stringsAsFactors = F)
  classes <- c("N")
  
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
  
  #
  # Map samples
  #
  
  mapSample <- function(sample){
    cd_doc()
    
    facets_data <- as.data.frame(read.table(paste("CSHL/Project_TUV_12995_B01_SOM_Targeted.2018-03-02/Sample_", sample, "/analysis/structural_variants/", sample, "--NA12878.cnv.facets.v0.5.2.txt", sep = ""), header = TRUE, sep="\t", stringsAsFactors=FALSE, quote=""))
    facets_data <- facets_data[,c(1, 10, 11, 5)]
    names(facets_data) <- c("chr", "start", "end", "cnlr")
    if (length(facets_data[facets_data$chr == "X", ]$chr) > 0) facets_data[facets_data$chr == "X", ]$chr <- "23" #TODO: May not work if multiple segments in X chromosome
    if (length(facets_data[facets_data$chr == "Y", ]$chr) > 0)facets_data[facets_data$chr == "Y", ]$chr <- "24"
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
    
    print(paste("Mapped ", sample, sep = ""))
    return(facets_data)
    
  }
  
  outputMappedSample <- function(sample, mapped_sample){
    cd_local()
    file_name <- paste("mappedFacetsFiles/", sample, "--NA12878.mapped.cnv.facets.v0.5.2.bed", sep = "")
    write.csv(mapped_sample, file = file_name, row.names = FALSE)
    print(paste("Wrote ", sample, " to ", file_name))
  }
  
  for(sample in loaded_samples){
    mapped_sample <- mapSample(sample)
    outputMappedSample(sample, mapped_sample)
  }
}

#
# For a dataframe from a single sample, map bins to bases
# - Knowing sample is potentially important since binning scheme may different
# - ^ This may be a challenge down the line since sample may be unknown in a concatenated datafram
# TODO: Do this once it is time for visualization
#
mapBinsToBases <- function() {
  
  
  
}

mapBasesToBins()
