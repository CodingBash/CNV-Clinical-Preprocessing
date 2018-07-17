setwd("~/Git-Projects/Git-Research-Projects/drug-response-prediction")
source("helperFunctions.R")
source("cnProfileVisualizationLibrary.R")



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
  xxFiles <- getAllXXFilesForTumor(tumorId)
  
  apply_func <- function(index){
    col_vec <- sapply(xxFiles, function(x) x[index, 4])
    return(sd(col_vec))
  }
  results <- lapply(seq(1, nrow(xxFiles[[tumorId]])), apply_func)
  result_df <- cbind(xxFiles[[tumorId]])
  result_df$value <- results
  
  # TODO: Next, organize code and apply MSE solution.
}

visualizeCNProfile(facets_snp_data = result_df, save = FALSE, ymin = 0, ymax = 1.5)
