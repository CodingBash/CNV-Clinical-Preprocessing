setwd("~/Git-Projects/Git-Research-Projects/drug-response-prediction")
source("helperFunctions.R")
source("cnProfileVisualizationLibrary.R")
library(BSgenome.Hsapiens.UCSC.hg19)
library(dplyr)

# TODO: May just pass back file
retrieveFacetsSegmentsFromObject <- function(tumorId, normalId, fitPrefix = "facetsG5Fit_", dir = "output/"){
  fit <- readRDS(paste(dir, fitPrefix, tumorId, "_", normalId, ".rds", sep = ""))
  return(fit$cncf)
}

segmentsToBedFormat <- function(facets_segment_data, median.clust = FALSE){
  facets_segment_data <- facets_segment_data[,c(1, 10, 11, if(median.clust == FALSE) 5 else 8)]
  facets_segment_data[[1]] <- paste("chr", facets_segment_data[[1]], sep="")
  names(facets_segment_data) <- c("chrom", "start", "end", "value")
  return(facets_segment_data)
}

for(tumorId in seq(1, 7)){
  all_tumor_segments <- data.frame()
  for(normalId in seq(1, 7)){
    if(normalId == tumorId) next
    try({
    cd_facets()
    raw_segments <- retrieveFacetsSegmentsFromObject(tumorId, normalId, fitPrefix = "facetsG5Fit_", dir = "output/prev_run_1/")
    bed <- segmentsToBedFormat(raw_segments)
    all_tumor_segments <- rbind(all_tumor_segments, bed)
    }, silent = TRUE)
  }
  visualizeCNProfile(all_tumor_segments,  save = FALSE)
}

getAllXXFilesForTumor <- function(tumorId){
  xxList <- list()
  cd_facets()
  for(normal.id in seq(1,7)){
    xx <- getFacetsXX(tumorId, normal.id, xxPrefix = "facetsG5XX_", dir = "output/prev_run_1/")
    snps <- snpsToBedFormat(xx$jointseg)
    xxList[[normal.id]] <- snps
  }
  return(xxList)
}

getAllFitFilesForTumor <- function(tumorId){
  fitList <- list()
  cd_facets()
  for(normal.id in seq(1,7)){
    fit <- getFacetsFit(tumorId, normal.id, fitPrefix = "facetsG5Fit_", dir = "output/prev_run_1/")
    segments <- segmentsToBedFormat(fit$cncf)
    fitList[[normal.id]] <- segments
  }
  return(fitList)
}

# Calculate SD
# TODO: Completed VERY fast, relook at this
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

# Calculate SD
# TODO: Completed VERY fast, relook at this
# - Not good to go SNP based due to inherent measurement variation. scratching
for(tumorId in seq(1, 7)){
  cd_facets()
  fitFiles <- getAllFitFilesForTumor(tumorId)
  
  nonmatching_segments <- data.frame()
  for(normalId in seq(1, 7)){
    if(tumorId == normalId){
      next
    }
    nonmatching_segments <- rbind(nonmatching_segments, fitFiles[[normalId]])
  }
  
  sizes <- generateChromosomeSizes(BSgenome.Hsapiens.UCSC.hg19)
  sizes <- data.frame(lapply(sizes, as.character), stringsAsFactors=FALSE)
  
  sizes[sizes$chrom == "chrX", ]$chrom <- "chr23"
  sizes[sizes$chrom == "chrY", ]$chrom <- "chr24"
  positions <- data.frame()
  for(sizes.i in seq(1, nrow(sizes))){
    positions <- rbind(positions, data.frame(chrom = sizes[sizes.i, 1], position = seq(1, sizes[sizes.i, 2], by=2500000)))
  }
  
  calculatedVariation <- function(index){
    return(sd(nonmatching_segments[nonmatching_segments$chrom == positions[index, ]$chrom & nonmatching_segments$start <= positions[index, ]$position &  nonmatching_segments$end >= positions[index, ]$position, 4]))
  }
  
  result <- lapply(seq(1, nrow(positions)), calculatedVariation)
  result_df <- cbind(positions)
  result_df$value <- result
  result_bed <- result_df[,c(1,2,2,3)]
  values <- as.vector(result_bed[[4]])
  values <- values[!is.na(values), ]
  # TODO: Next, organize code and apply MSE solution.
  visualizeCNProfile(facets_segment_data = fitFiles[[tumorId]],facets_snp_data = result_bed, save = FALSE, ymin = 0, ymax = 1.5)
  # TODO: get mean, visualize better, add segments to visualization, MSE
}

visualizeCNProfile(facets_snp_data = result_bed, save = FALSE, ymin = 0, ymax = 1.5)
