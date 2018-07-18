#
# This script contains a set of helper functions for the other CNA scripts
#

#
# Sets the directory to the project workspace
#
cd_local <- function() {
  setwd("~/Git-Projects/Git-Research-Projects/drug-response-prediction") # TODO: cd not working in notebook - home directory different
}

#
# Sets the directory to facets_testing workspace
# TODO: all scripts that go under facets_testing should be moved to that workspace
#
cd_facets <- function(subdir = ""){
  setwd(paste0("~/Git-Projects/Git-Research-Projects/FACETS_nonmatching_test/", subdir))
}


#
# Sets the directory to the project workspace
#
cd_doc <- function() {
  setwd("~") # TODO: cd not working in notebook - home directory different
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
# Retrieve a sample's FACETS segmented cluster (cnlr.median.clust) data from specified directory from CNprep results
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
segmentsToBedFormat <- function(facets_segment_data, median.clust = FALSE){
  facets_segment_data <- facets_segment_data[,c(1, 10, 11, if(median.clust == FALSE) 5 else 8)]
  facets_segment_data[[1]] <- paste("chr", facets_segment_data[[1]], sep="")
  names(facets_segment_data) <- c("chrom", "start", "end", "value")
  return(facets_segment_data)
}

#
# Simplifies FACETS segment original format into BED format with
# the columns: "chrom", "chrom.start", "chrom.end", "median cnlr"
#
segmentsToBedFormatWithClusters <- function(facets_segment_data){
  facets_segment_data <- facets_segment_data[,c(1, 10, 11, 8)]
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
  facets_snp_data[[1]] <- paste("chr", facets_snp_data[[1]], sep="")
  names(facets_snp_data) <- c("chrom", "start", "end", "value")
  return(facets_snp_data)
}

#
# Concatenates all segments in sample that follow a specific set of events 
#
# @param event - either "A" (amplification) or "D" (deletion) segments or "N" (neutral) segments are extracted
# @param samples - names of samples to retrieve
# @param chromosomeSizes - df of chromosomeSizes for rescaling (if necessary)
# @param dir - directory of the input files
# @param extension - extension of the input files
# TODO: segment filename still too hardcoded. Allow caller to send file name themselves
# @param rescaleInput - indicator if sample segments should be scaled (usually if it is not absolute scaled and is chromsome scaled instead). If TRUE, rescaling with chromosomeSizes is performed
# @param inSampleFolder - ad-hoc solution when the sample file is contained in a folder with sample name
# @param ampCall - lower threshold to call amplifications
# @param delCall - higher threshold to call deletions
#
# TODO: This method is also used in segmentClustering.R script. Perhaps move this function to a more general library?
# 
# @DEPRECTATED use selectSegmentsWithEvents, which modularizes the loading of segments and the subsetting of all segments
selectSegmentsWithEventsDEPRECATED <- function(events, samples, chromosomeSizes, dir, extension = "cnv.facets.v0.5.2.txt", inSampleFolder = FALSE, rescaleInput = FALSE, ampCall = 0.2, delCall = -0.235){
  ###############################################
  # DEPRECATED - user selectSegmentsWithEvents  #
  ###############################################
  
  totalSelectedSegments <- data.frame()
  
  loaded_segments <- list(NA)
  loaded_segments.index <- 1
  for(sample in samples){
    segments <- as.data.frame(read.table(paste(dir, if(inSampleFolder == TRUE) paste("Sample_", sample, "/analysis/structural_variants/", sep = "") else "",sample, "--NA12878.", extension, sep = ""), header = TRUE, sep="\t", stringsAsFactors=FALSE, quote=""))
    selected_segments <- data.frame()
    
    # Retrieve amplification events from sample if asked
    if("A" %in% events){
      selected_segments <- rbind(selected_segments, segments[segments$X.cnlr.median. > ampCall,])  
    }
    
    # Retrieve deletion events from sample if asked
    if ("D" %in% events){
      selected_segments <- rbind(selected_segments, segments[segments$X.cnlr.median. < delCall,])
    }
    
    # Retrieve nuetral events from sample if asked
    if("N" %in% events){
      selected_segments <- rbind(selected_segments, segments[segments$X.cnlr.median. >= delCall & segments$X.cnlr.median. <= ampCall,]) # TODO: This is an untested line of code
    }
    
    # If any segments from sample selected, let's preprocess the dataframe and add to total list
    if(nrow(selected_segments) != 0){
      # Filters to only chrom, start, end, cnlr
      selected_segments <- selected_segments[,c(1, 10, 11, 5)] # TODO: Just added CNLR, may have issues in CORE script
      
      names(selected_segments) <- c("chrom", "start", "end", "cnlr")
      
      if(rescaleInput == TRUE){
        selected_segments <- chromsomeToAbsoluteBPConversion(selected_segments, chromosomeSizes)
      }
      
      totalSelectedSegments <- rbind(totalSelectedSegments, selected_segments)
    }
  }
  
  returnme <-  cbind(totalSelectedSegments)
  returnme <- returnme[returnme$chrom != "X" & returnme$chrom != "Y",] # REMOVE X and Y chromosome
  returnme$chrom <- as.numeric(returnme$chrom)
  return(returnme)
}

#
# Concatenates all segments in sample that follow a specific set of events 
#
# @param event - either "A" (amplification) or "D" (deletion) segments or "N" (neutral) segments are extracted
# @param samples - names of samples to retrieve
# @param chromosomeSizes - df of chromosomeSizes for rescaling (if necessary)
# @param dir - directory of the input files
# @param extension - extension of the input files
# TODO: segment filename still too hardcoded. Allow caller to send file name themselves
# @param rescaleInput - indicator if sample segments should be scaled (usually if it is not absolute scaled and is chromsome scaled instead). If TRUE, rescaling with chromosomeSizes is performed
# @param inSampleFolder - ad-hoc solution when the sample file is contained in a folder with sample name
# @param ampCall - lower threshold to call amplifications
# @param delCall - higher threshold to call deletions
#
# TODO: This method is also used in segmentClustering.R script. Perhaps move this function to a more general library?
#
selectSegmentsWithEvents <- function(events, samples, chromosomeSizes, dir, extension = "cnv.facets.v0.5.2.txt", inSampleFolder = FALSE, rescaleInput = FALSE, ampCall = 0.2, delCall = -0.235){
  # TODO: May have trouble supporting missing or default parameters
  segmentList <- retrieveSegmentListFromSamples(samples, dir, extension, inSampleFolder)
  selectedSegments <- subsetAllSegmentsByEvent(segmentList, events, chromosomeSizes, rescaleInput, ampCall, delCall)
  return(selectedSegments)
}

#
# From a list of samples, retrieve the segments in a key-value list (where K is the sample name, and V is the segment dataframe)
#
retrieveSegmentListFromSamples <- function(samples, dir, extension = "cnv.facets.v0.5.2.txt", inSampleFolder = FALSE){
  segmentList <- list(NA)
  for(sample in samples){
    segments <- as.data.frame(read.table(paste(dir, if(inSampleFolder == TRUE) paste("Sample_", sample, "/analysis/structural_variants/", sep = "") else "",sample, "--NA12878.", extension, sep = ""), header = TRUE, sep="\t", stringsAsFactors=FALSE, quote=""))
    segmentList[[sample]] <- segments
  }
  return(segmentList)
}

#
# From a list of segments, subset the segments to only include amplfication, deletion, and/or neutral events based on user specifications
#
subsetAllSegmentsByEvent <- function(segmentList, events, chromosomeSizes, rescaleInput = FALSE, ampCall = 0.2, delCall = -0.235){
  totalSelectedSegments <- data.frame()
  names(segmentList) <- as.character(seq_along(segmentList)) # Convert indices to character names so we can consistently iterate through any list by name
  for(segmentName in names(segmentList)){
    segments <- segmentList[[segmentName]]
    if(is.null(segments)){
      next
    }
    selected_segments <- data.frame()
    # Retrieve amplification events from sample if asked
    if("A" %in% events){
      selected_segments <- rbind(selected_segments, segments[segments[[5]] > ampCall,])  
    }
    
    # Retrieve deletion events from sample if asked
    if ("D" %in% events){
      selected_segments <- rbind(selected_segments, segments[segments[[5]] < delCall,])
    }
    
    # Retrieve nuetral events from sample if asked
    if("N" %in% events){
      selected_segments <- rbind(selected_segments, segments[segments[[5]] >= delCall & segments[[5]] <= ampCall,]) # TODO: This is an untested line of code
    }
    
    # If any segments from sample selected, let's preprocess the dataframe and add to total list
    if(nrow(selected_segments) != 0){
      # Filters to only chrom, start, end, cnlr
      selected_segments <- selected_segments[,c(1, 10, 11, 5)] # TODO: Just added CNLR, may have issues in CORE script
      
      names(selected_segments) <- c("chrom", "start", "end", "cnlr")
      
      if(rescaleInput == TRUE){
        selected_segments <- chromsomeToAbsoluteBPConversion(selected_segments, chromosomeSizes)
      }
      
      totalSelectedSegments <- rbind(totalSelectedSegments, selected_segments)
    }
  }
  
  returnme <-  cbind(totalSelectedSegments)
  returnme <- returnme[returnme$chrom != "X" & returnme$chrom != "Y",] # REMOVE X and Y chromosome
  returnme$chrom <- as.numeric(returnme$chrom)
  return(returnme)
}
#
# Given a genome (i.e. hg19), generate the chromosome sizes
#
generateChromosomeSizes <- function(genome){
  seqlengths(genome) <- seqlengths(Hsapiens)
  
  # Create chromosome vector
  chrom_vec <- c(NA)
  chrom_vec.index <- 1
  for(i in append(seq(1,22, by=1), c("X", "Y"))){
    chrom_vec[chrom_vec.index] <- paste("chr", i, sep = "")  
    chrom_vec.index <- chrom_vec.index + 1
  }
  chromosomeSizes <- data.frame(stringsAsFactors = FALSE)
  
  for(chrom_i in chrom_vec){
    df = data.frame(chrom = chrom_i, size = seqlengths(genome)[chrom_i])
    chromosomeSizes <- rbind(chromosomeSizes, df)
  }
  return(chromosomeSizes)
}

chromsomeToAbsoluteBPConversionForSingleEntry <- function(chrom, start, end, chromosomeSizes){
  chrom_r <- chrom
  if (is.na(chrom_r) || length(chrom_r) == 0){
    next # TODO: This is to resolve the NA row. Where did it come from?
  } 
  total_bp <- 0
  if(chrom_r %in% seq(2,22)){
    for(i in seq(1, as.numeric(chrom_r) - 1)){
      total_bp <- total_bp + chromosomeSizes[paste("chr", i, sep = ""), ]$size
    }  
  } else if (chrom_r == "X") {
    for(i in seq(1, 22)){
      total_bp <- total_bp + chromosomeSizes[paste("chr", i, sep = ""), ]$size
    }  
  }  else if (chrom_r == "Y") {
    for(i in seq(1, 22)){
      total_bp <- total_bp + chromosomeSizes[paste("chr", i, sep = ""), ]$size
    }  
    total_bp <- total_bp + chromosomeSizes["chrX", ]$size
  }
  abs_start <- start + total_bp
  abs_end <- end + total_bp
  returnme <- data.frame(start = abs_start, end = abs_end)
  return(returnme)
}

#
# Take a input of multiple segments with the chromosomeSizes, and convert
# segment maploc from chrom.location to absolute.location in bp units
#
chromsomeToAbsoluteBPConversion <- function(input, chromosomeSizes){
  for(row.index in seq(1, nrow(input))){
    # Rescale row
    absoluteRow <- chromsomeToAbsoluteBPConversionForSingleEntry(chrom = input[row.index, ]$chrom, start = input[row.index, ]$start, end = input[row.index, ]$end, chromosomeSizes = chromosomeSizes)
    
    #
    # Update start and end maploc with new absolute location
    #
    input[row.index, ]$start <- absoluteRow$start
    input[row.index, ]$end <- absoluteRow$end
  }
  return(input)
}

#
# From a dataframe, create a TSV file (bed file)
#
createBedFile <- function(objectToWrite, filename){
  write.table(objectToWrite, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE, file = filename)
}