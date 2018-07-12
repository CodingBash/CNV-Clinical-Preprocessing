#
# This script contains functions to display copy number profiles
#


#
# Import libraries
# WARNING: Gtrellis must be installed from bioconductor before use
#
library(gtrellis)
library(circlize)
library(ComplexHeatmap)

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


#
# Visualizes CN profile from segment/SNP data using GTRELLIS visualization library
#
visualizeCNProfile <- function(facets_segment_data, facets_snp_data, categories, save = FALSE, saveDir = "", saveMeta = ""){
  if(save == TRUE){
    pdf(paste(saveDir, "cnprofile_", saveMeta, "_", sample, ".pdf", sep = ""), width = 44, height = 16)
  }
  # Add legend
  lgd <- Legend(at = c("duplication", "nuetral", "deletion"), title = "Class", type = "lines", legend_gp = gpar(col = c("orange", "blue", "red")))
  
  #
  # Set GTRELLIS layout
  #
  track_height <- c(1,8,1)
  track_axis <- c(FALSE, TRUE, FALSE)
  track_ylim <- range(seq(-6,2))
  nrow <- 3
  n_track <- 3
  byrow <- FALSE
  species <- "hg19"
  legend <- lgd
  if(missing(categories)){
    gtrellis_layout(track_height = track_height,
                    track_axis = track_axis, 
                    track_ylim = track_ylim, 
                    nrow = nrow, 
                    n_track = n_track, 
                    byrow = byrow, 
                    species = species, 
                    legend = legend)
  } else {
    gtrellis_layout(track_height = track_height,
                    track_axis = track_axis, 
                    track_ylim = track_ylim, 
                    nrow = nrow, 
                    n_track = n_track, 
                    byrow = byrow, 
                    species = species, 
                    category = categories,
                    legend = legend)
  }
  
  # Add top chromosome labels
  add_track(panel_fun = function(gr) {
    chr = get_cell_meta_data("name")  
    grid.rect(gp = gpar(fill = "#EEEEEE"))
    grid.text(chr)
  })
  
  # Add SNP/bin track to GTRELLIS  
  add_points_track(facets_snp_data, facets_snp_data$value, gp = gpar(fill = "gray"))
  # Add segment track overlayed over the SNP/bins to GTRELLIS
  add_segments_track(facets_segment_data, facets_segment_data$value, track = current_track(), gp = gpar(col = ifelse(facets_segment_data$value > 0.2, "orange", ifelse(facets_segment_data$value > -0.23, "blue", "red")), lwd = 4))
  
  # Add cytobands to bottom of panel to GTRELLIS
  cytoband_df = circlize::read.cytoband(species = "hg19")$df
  add_track(cytoband_df, panel_fun = function(gr) {
    cytoband_chr = gr
    grid.rect(cytoband_chr[[2]], unit(0, "npc"),
              width = cytoband_chr[[3]] - cytoband_chr[[2]], height = unit(1, "npc"),
              default.units = "native", hjust = 0, vjust = 0,
              gp = gpar(fill = circlize::cytoband.col(cytoband_chr[[5]])))
    grid.rect(min(cytoband_chr[[2]]), unit(0, "npc"),
              width = max(cytoband_chr[[3]]) - min(cytoband_chr[[2]]), height = unit(1, "npc"),
              default.units = "native", hjust = 0, vjust = 0,
              gp = gpar(fill = "transparent"))
  })
  
  if(save == TRUE){
    dev.off()
  }
}
