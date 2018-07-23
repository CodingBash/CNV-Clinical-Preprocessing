#
# Visualize CN profile that is color coded based on cluster assignment from CNprep
#

#
# This script produces the CN profile images and saves in preferred directory
#

#
# Load sources
# WARNING: If sources don't load, manually setwd() to src file location
#
source("cnProfileVisualizationLibrary.R")
source("helperFunctions.R")
source("visualizeGMMLibrary.R")

generateBaseAxis <- function(chromosomeSizes){
  base_axis <- data.frame()
  for(chromosomeSizes.index in seq_along(chromosomeSizes$size)){
    entry <- data.frame(chrom = chromosomeSizes[chromosomeSizes.index, ]$chrom, start = 1, end = chromosomeSizes[chromosomeSizes.index, ]$size, value = 0.0)
    base_axis <- rbind(base_axis, entry)
  }
  return(base_axis)
}


visualizeClusterAssignment <- function(sample, base_axis, categories, segDir = "segClusteringResults/prev_run_7_19_2018_2/", save = FALSE, saveDir = "./clusterAssignmentVisualization/"){
  #
  # Retrieve and format input data
  #
  facets_segment_clusters <- retrieveSegtable(sample, segDir)
  gaussian_comps <- unique(facets_segment_clusters[,c("maxzmean", "maxzsigma")])
  row.names(gaussian_comps) <- seq_along(row.names(gaussian_comps))
  component_indices <- unlist(lapply(seq(1, nrow(facets_segment_clusters)), function(index){
    # Return the index of the GMM component that the segment belongs to
    return(as.numeric(rownames(gaussian_comps[gaussian_comps$maxzmean == facets_segment_clusters[index, ]$maxzmean & gaussian_comps$maxzsigma == facets_segment_clusters[index, ]$maxzsigma, ])))
  }))

  # Set the segment index as the corrensponding name of each component index
  facets_segment_clusters_bed <- segmentClustersToBedFormat(facets_segment_clusters, value = 5)
  indices <- unlist(lapply(seq(1, nrow(facets_segment_clusters)), function(index){
    if(facets_segment_clusters[index, ]$chrom.pos.start %in% facets_segment_clusters_bed$start & facets_segment_clusters[index, ]$chrom.pos.end %in% facets_segment_clusters_bed$end){
      return(index)
    }
  }))
  
  
  component_indices <- component_indices[indices]
  names(component_indices) <- rownames(facets_segment_clusters_bed)

  if(!missing(base_axis)){
    prev_length <- nrow(facets_segment_clusters_bed)
    facets_segment_clusters_bed <- rbind(facets_segment_clusters_bed, base_axis)
    axis_indices <- rep(6, nrow(base_axis))
    names(axis_indices) <- seq(prev_length+ 1, nrow(facets_segment_clusters_bed))
    component_indices <- append(component_indices, axis_indices)
  }
  
  #
  # Generate CN profile from input data
  #
  visualizeCNProfile(facets_segment_data = facets_segment_clusters_bed, categories = categories, color_id = component_indices, save = save, saveDir = saveDir, saveMeta = paste0("clusterassignment_", sample), title=paste0("CNprep Cluster Assignment for sample=", sample))  
}

