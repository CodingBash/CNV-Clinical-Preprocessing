setwd("~/Git-Projects/Git-Research-Projects/drug-response-prediction")
source("helperFunctions.R")
source("cnProfileVisualizationLibrary.R")
source("facetsAnalysisLibrary.R")
library(BSgenome.Hsapiens.UCSC.hg19)


#
# Visualize all segments per each tumor overlapped among eachother
#
for(tumorId in seq(1, 7)){
  all_tumor_segments <- data.frame()
  for(normalId in seq(1, 7)){
    #if(normalId == tumorId) next
    try({
    cd_facets()
    raw_segments <- retrieveFacetsSegmentsFromObject(tumorId, normalId, fitPrefix = "facetsG5Fit_", dir = "output/prev_run_1/")
    bed <- segmentsToBedFormat(raw_segments)
    all_tumor_segments <- rbind(all_tumor_segments, bed)
    }, silent = TRUE)
  }
  visualizeCNProfile(all_tumor_segments,  save = FALSE)
}


# Calculate SD across each SNP in the nonmatching samples
# TODO: Completed VERY fast, relook at this
# @deprecated - very inefficient code to calculate SD
#
for(tumorId in seq(1, 7)){
  cd_facets()
  xxFiles <- getAllXXFilesForTumor(tumorId)
  bed <- xxFiles[[tumorId]]
  sdBed <- data.frame()
  for(bed.i in seq(1, nrow(bed))){
    bed.i.sum <- 0
    for(normalId in seq(1, 7)){
      if(normalId == tumorId) next
      bed.i.sum <- xxFiles[[normalId]][bed.i,]$value
    }
    mean <- bed.i.sum / 6.0 # TODO: 6 (sample size) should not be hardcoded
    sum_sq_diff <- 0
    for(normalId in seq(1, 7)){
      if(normalId == tumorId) next
      sum_sq_diff <- (xxFiles[[normalId]][bed.i,]$value - mean) ^ 2
    }
    sd <- sum_sq_diff / 6.0
    entry <- data.frame(chrom = bed[bed.i, ]$chrom, start = bed[bed.i,]$start, end = bed[bed.i, ]$end, value = mse)
    sdBed <- rbind(sdBed, entry)
  }    
}

# Calculate SD across each SNP in the nonmatching samples
# @deprecated Not good to go SNP based due to inherent measurement variation. Scratching
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
  values <- as.vector(result_bed[[4]])
  values <- values[!is.na(values), ]
  
  # Visualize results
  visualizeCNProfile(facets_segment_data = fitFiles[[tumorId]],facets_snp_data = result_bed, save = FALSE, ymin = 0, ymax = 1.5)
  
  # TODO: Next, organize code and apply MSE solution.
  # TODO: get mean, improve visualization, add segments to visualization, MSE
}



