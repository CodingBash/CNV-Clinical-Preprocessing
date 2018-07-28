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
# TODO: Have documentation on each parameter
# TODO: Allow user to send in legend info (legend is hardcoded as of now)
#
visualizeCNProfile <- function(facets_segment_data, facets_snp_data, line_data, Acores, Dcores, categories, save = FALSE, saveDir = "", saveMeta = "", ymin = -6, ymax = 2, title = "", color_id){
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
    gtrellis_layout(title = title,
                    track_height = track_height,
                    track_axis = track_axis, 
                    track_ylim = track_ylim, 
                    nrow = nrow, 
                    n_track = n_track, 
                    byrow = byrow, 
                    species = species, 
                    legend = legend)
  } else {
    gtrellis_layout(title = title,
                    track_height = track_height,
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
  add_track(panel_fun = function(color_id) {
    chr = get_cell_meta_data("name")  
    grid.rect(gp = gpar(fill = "#EEEEEE"))
    grid.text(chr)
  })
  
  colors <- c("blue", "red", "orange", "green", "purple", "black", "gray")
  determineSegmentColors <- function(facets_segment_data, color_id){
    # TODO: Build determineColor functions
    if(missing(color_id)){
      col = ifelse(facets_segment_data$value > 0.2, "orange", ifelse(facets_segment_data$value > -0.23, "blue", "red"))
    } else {
      col = unlist(lapply(seq_along(rownames(facets_segment_data)), function(index){
        color_val <- color_id[[rownames(facets_segment_data)[[index]]]] 
        return(colors[[color_val]])
      }))
    }
    return(col)
  }
  
  #
  # TODO: The logic in the following conditions may be incorrect.
  # Below, I am attempting to overlap 3 different pieces of data (segment, SNP, lines)
  # depending on what is and is not provided (since I need to use track = current_track()) if
  # the track is already preoccupied by one of the other 3 pieces of data
  #
  if(!missing(facets_segment_data)){
    #add_segments_track(facets_segment_data, facets_segment_data$value, gp = gpar(col = ifelse(facets_segment_data$value > 0.2, "orange", ifelse(facets_segment_data$value > -0.23, "blue", "red")), lwd = 2))
    if(!missing(facets_snp_data)){
      add_points_track(facets_snp_data, facets_snp_data$value, gp = gpar(fill = "gray"))  
      if(!missing(line_data)){
        add_lines_track(line_data, line_data[[4]], gp = gpar(col = "purple", lwd = 5), track = current_track())
      }
      add_segments_track(facets_segment_data, facets_segment_data$value, gp = gpar(col = determineSegmentColors(facets_segment_data, color_id), lwd = 2), track = current_track())
    } else {
      if(!missing(line_data)){
        add_lines_track(line_data, line_data[[4]], gp = gpar(col = "purple", lwd = 5), track = current_track())
        add_segments_track(facets_segment_data, facets_segment_data$value, gp = gpar(col = determineSegmentColors(facets_segment_data, color_id), lwd = 2), track = current_track())
      }
      add_segments_track(facets_segment_data, facets_segment_data$value, gp = gpar(col = determineSegmentColors(facets_segment_data, color_id), lwd = 2))
    }
  } else if(!missing(facets_snp_data)){
    add_points_track(facets_snp_data, facets_snp_data$value, gp = gpar(fill = "gray"))  
    if(!missing(line_data)){
      add_lines_track(line_data, line_data[[4]], gp = gpar(col = "purple", lwd = 5), track = current_track())
    }
  } else if(!missing(line_data)){
    add_lines_track(line_data, line_data[[4]], gp = gpar(col = "purple", lwd = 5))
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


