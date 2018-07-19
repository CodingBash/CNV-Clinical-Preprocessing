#
# This script calculates the SD and RMSE of the TCGA FACETS output from the 
# nonmatching tumor-normal pairs. The statistics were calculate per SNP and
# by segments (sampled at intervals)
#

#
# Import source libraries
#
source("helperFunctions.R")
source("facetsAnalysisLibrary.R")
library(BSgenome.Hsapiens.UCSC.hg19)

# TODO: Code can be modularized

#
# Calculate SD across each SNP in the nonmatching samples
#
calculateSnpSD <- function(tumorId){
  xxFiles <- getAllXXFilesForTumor(tumorId)
  bed <- xxFiles[[tumorId]]
  calculateSnpVariation <- function(index){
    snp_logR <- unlist(lapply(seq(1, 7), function(normalId){
      if(tumorId != normalId){
        val <- xxFiles[[normalId]][index,]$value
        if(!is.na(val)) return(val)
      }
    }))
    return(sd(snp_logR))
  }
  
  result <- unlist(lapply(seq(1, nrow(bed)), calculateSnpVariation))
  mean_sd <- mean(result[!is.na(result)])
  print(paste0("The mean SD for SNPs of tumorID=", tumorId, " is SD=", mean_sd))
  
  sd_bed <- cbind(bed)
  sd_bed$value <- result
  
  return(sd_bed)
}


#
# Calculate RMSE across each SNP in the nonmatching samples
#
calculateSnpRMSE <- function(tumorId) {
  xxFiles <- getAllXXFilesForTumor(tumorId)
  bed <- xxFiles[[tumorId]]
  calculateSnpRMSE <- function(index){
    squared_errors <- unlist(lapply(seq(1, 7), function(normalId){
      if(tumorId != normalId){
        squared_error <- (xxFiles[[normalId]][index,]$value - xxFiles[[tumorId]][index,]$value) ^ 2
        if(!is.na(squared_error)) return(squared_error)
      }
    }))
    return(sqrt(mean(squared_errors)))
  }
  
  result <- unlist(lapply(seq(1, nrow(bed)), calculateSnpRMSE))
  
  mean_rmse <- mean(result[!is.na(result)])
  print(paste0("The mean RMSE for SNPs of tumorId=", tumorId, " is MRMSE=", mean_rmse))
  
  rmse_bed <- cbind(bed)
  rmse_bed$value <- result
  return(rmse_bed)
}


#
# Calculate SD across each segment in the nonmatching samples
#
calculateSegmentSD <- function(tumorId){
  # Get fit files as BED format
  fitFiles <- getAllFitFilesForTumor(tumorId)
  
  # Get nonmatching fit segments
  nonmatching_segments <- data.frame()
  for(normalId in seq(1, 7)){
    if(tumorId == normalId){
      next
    }
    nonmatching_segments <- rbind(nonmatching_segments, fitFiles[[normalId]])
  }
  
  # TODO: This should only be calculated once (when iterating through tumorIds)
  # Get sizes to know what positions to sample variation at
  sizes <- generateChromosomeSizes(BSgenome.Hsapiens.UCSC.hg19)
  sizes <- data.frame(lapply(sizes, as.character), stringsAsFactors=FALSE)
  sizes[sizes$chrom == "chrX", ]$chrom <- "chr23"
  sizes[sizes$chrom == "chrY", ]$chrom <- "chr24"
  positions <- data.frame()
  for(sizes.i in seq(1, nrow(sizes))){
    # Sample position at every n interval
    positions <- rbind(positions, data.frame(chrom = sizes[sizes.i, 1], position = seq(1, sizes[sizes.i, 2], by=2500000))) 
  }
  
  # Calculate standard deviation at each position
  calculatedVariation <- function(index){
    return(sd(nonmatching_segments[nonmatching_segments$chrom == positions[index, ]$chrom & nonmatching_segments$start <= positions[index, ]$position &  nonmatching_segments$end >= positions[index, ]$position, 4]))
  }
  result <- unlist(lapply(seq(1, nrow(positions)), calculatedVariation))
  mean_sd <- mean(result[!is.na(result)])
  print(paste0("The mean SD for segments of tumorID=", tumorId, " is SD=", mean_sd))
  
  
  #
  # Format results
  #
  result_df <- cbind(positions)
  result_df$value <- result
  sd_bed <- result_df[,c(1,2,2,3)]
  return(sd_bed)
}

#
# Calculate RMSE across each segment in the nonmatching samples
#
calculateSegmentRMSE <- function(tumorId) {
  # Get fit files as BED format
  fitFiles <- getAllFitFilesForTumor(tumorId)
  
  # Get nonmatching fit segments
  nonmatching_segments <- data.frame()
  for(normalId in seq(1, 7)){
    if(tumorId == normalId){
      next
    }
    nonmatching_segments <- rbind(nonmatching_segments, fitFiles[[normalId]])
  }
  
  # TODO: This should only be calculated once (when iterating through tumorIds)
  # Get sizes to know what positions to sample variation at
  sizes <- generateChromosomeSizes(BSgenome.Hsapiens.UCSC.hg19)
  sizes <- data.frame(lapply(sizes, as.character), stringsAsFactors=FALSE)
  sizes[sizes$chrom == "chrX", ]$chrom <- "chr23"
  sizes[sizes$chrom == "chrY", ]$chrom <- "chr24"
  positions <- data.frame()
  for(sizes.i in seq(1, nrow(sizes))){
    # Sample position at every n interval
    positions <- rbind(positions, data.frame(chrom = sizes[sizes.i, 1], position = seq(1, sizes[sizes.i, 2], by=2500000))) 
  }
  
  matching_segment <- fitFiles[[tumorId]]
  
  # Calculate standard deviation at each position
  calculateRMSE <- function(index){
    
    overlapped_nm_segments <- nonmatching_segments[nonmatching_segments$chrom == positions[index, ]$chrom & nonmatching_segments$start <= positions[index, ]$position &  nonmatching_segments$end >= positions[index, ]$position, 4]
    overlapped_m_segment <- matching_segment[matching_segment$chrom == positions[index, ]$chrom & matching_segment$start <= positions[index, ]$position &  matching_segment$end >= positions[index, ]$position, 4]
    
    squared_errors <- unlist(lapply(seq_along(overlapped_nm_segments), function(nm_index){
      squared_error <- (overlapped_nm_segments[[nm_index]] - overlapped_m_segment) ^ 2
    }))
    
    rmse <- sqrt(mean(squared_errors)) 
    return(rmse)
  }
  result <- unlist(lapply(seq(1, nrow(positions)), calculateRMSE))
  
  mean_rmse <- mean(result[!is.na(result)])
  print(paste0("The mean RMSE for segments of tumorId=", tumorId, " is MRMSE=", mean_rmse))
  
  #
  # Format results
  #
  result_df <- cbind(positions)
  result_df$value <- result
  rmse_bed <- result_df[,c(1,2,2,3)]
  return(rmse_bed)
}