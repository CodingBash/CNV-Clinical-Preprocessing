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
library(RColorBrewer)

#
# Visualizes CN profile from segment/SNP data using GTRELLIS visualization library
#
visualizeCNProfile <- function(facets_segment_data, facets_snp_data, Acores, Dcores, categories, save = FALSE, saveDir = "", saveMeta = "", ymin = -6, ymax = 2){
  if(save == TRUE){
    pdf(paste(saveDir, "cnprofile_", saveMeta, ".pdf", sep = ""), width = 44, height = 16)
  }
  # Add legend
  lgd <- Legend(at = c("duplication", "nuetral", "deletion"), title = "Class", type = "lines", legend_gp = gpar(col = c("orange", "blue", "red")))
  
  #
  # Set GTRELLIS layout
  #
  track_height <- if(!missing(Acores) || !missing(Dcores)) c(1,8,2,1) else c(1,8,1) 
  track_axis <-  if(!missing(Acores) || !missing(Dcores)) c(FALSE, TRUE, FALSE, FALSE) else c(FALSE, TRUE, FALSE)
  track_ylim <- range(seq(ymin, ymax))
  nrow <- 3
  n_track <- if(!missing(Acores) || !missing(Dcores)) 4 else 3
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
  
  if(!missing(facets_segment_data)){
    # Add segment track
    add_segments_track(facets_segment_data, facets_segment_data$value, gp = gpar(col = ifelse(facets_segment_data$value > 0.2, "orange", ifelse(facets_segment_data$value > -0.23, "blue", "red")), lwd = 2))
  }
  # If SNPs were provided
  if(!missing(facets_snp_data)){
    # Add SNP/bin track to GTRELLIS (overlayed over segment)
    if(!missing(facets_segment_data)){
      add_points_track(facets_snp_data, facets_snp_data$value, track = current_track(), gp = gpar(fill = "gray"))  
    } else {
      add_points_track(facets_snp_data, facets_snp_data$value, gp = gpar(fill = "gray"))  
    }
  }
  
  if(!missing(Acores) || !missing(Dcores)) {
    if(!missing(Acores) && !missing(Dcores)) {
      add_segments_track(Acores, 0, gp = gpar(col = "orange", lwd = 15))
      add_segments_track(Dcores, -2, gp = gpar(col = "red", lwd = 15), track = current_track())
    } else if(!missing(Acores)){
      add_segments_track(Acores, 0, gp = gpar(col = "orange", lwd = 15))
    } else {
      add_segments_track(Dcores, -2, gp = gpar(col = "red", lwd = 15))
    }
  }
  
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


