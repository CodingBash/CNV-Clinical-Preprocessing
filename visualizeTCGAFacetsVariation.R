setwd("~/Git-Projects/Git-Research-Projects/drug-response-prediction")
source("helperFunctions.R")
source("cnProfileVisualizationLibrary.R")
source("facetsAnalysisLibrary.R")
library(BSgenome.Hsapiens.UCSC.hg19)

# Calculate SD across each SNP in the nonmatching samples
#
for(tumorId in seq(1, 7)){
  cd_facets()
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
  
  sd_bed <- cbind(bed)
  sd_bed$value <- result
  
  visualizeCNProfile(facets_snp_data = sd_bed, categories = c("chr14"), save = FALSE, ymin = 0, ymax = 1.5)
}

# Calculate MSE across each SNP in the nonmatching samples
#
for(tumorId in seq(1, 7)){
  cd_facets()
  xxFiles <- getAllXXFilesForTumor(tumorId)
  bed <- xxFiles[[tumorId]]
  calculateSnpMSE <- function(index){
    squared_errors <- unlist(lapply(seq(1, 7), function(normalId){
      if(tumorId != normalId){
        squared_error <- (xxFiles[[normalId]][index,]$value - xxFiles[[tumorId]][index,]$value) ^ 2
        if(!is.na(squared_error)) return(squared_error)
      }
    }))
    return(mean(squared_errors))
  }
  
  result <- unlist(lapply(seq(1, nrow(bed)), calculateSnpMSE))
  
  mse_bed <- cbind(bed)
  mse_bed$value <- result
  
  visualizeCNProfile(facets_snp_data = mse_bed, categories = c("chr14"), save = FALSE, ymin = 0, ymax = 1.5)
}

#
# Calculate SD across each segment in the nonmatching samples
#
for(tumorId in seq(1, 7)){
  cd_facets()
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
  result <- lapply(seq(1, nrow(positions)), calculatedVariation)
  
  #
  # Format results
  #
  result_df <- cbind(positions)
  result_df$value <- result
  result_bed <- result_df[,c(1,2,2,3)]
  # Visualize results
  visualizeCNProfile(facets_snp_data = result_bed, categories = c("chr14"), save = FALSE, ymin = 0, ymax = 1)
  
  # TODO: Next, organize code and apply MSE solution.
  # TODO: get mean, improve visualization, add segments to visualization, MSE
}

