library(gtrellis)
library(circlize)
library(ComplexHeatmap)

# TODO: Verify outputCOREobj results, and create Jupyter notebook for CORE visualization

densityInput <- cbind(absoluteToChromosomeBPConversion(inputCORESegments, chromosomeSizes))
densityInput$chrom <- paste("chr", densityInput$chrom, sep = "")
dataInputCORE_density = circlize::genomicDensity(densityInput, window.size = 1e7, overlap = FALSE)

lgd = Legend(at = c("duplication", "nuetral", "deletion"), title = "Class", type = "lines", legend_gp = gpar(col = c("orange", "blue", "red")))
gtrellis_layout(track_height = unit.c(2*grobHeight(textGrob("chr1")), 
                                      unit(1.0, "null"), 
                                      unit(0.5, "null"), 
                                      unit(3, "mm")),
                track_axis = c(FALSE, TRUE, TRUE, FALSE), 
                track_ylim = c(0, 1, range(data.frame(COREtable$score)), c(0, max(dataInputCORE_density[[4]])), 0, 1),
                nrow = 3, 
                n_track = 4, 
                byrow = FALSE, 
                species="hg19",
                legend = lgd
                #,category = c("chr19")
)

# gtrellis_layout(track_height = c(2,5,1),
#                 track_axis = c(FALSE, TRUE, FALSE), 
#                 track_ylim = range(data.frame(facets_bins_data$X.cnlr.)),
#                 nrow = 3, 
#                 n_track = 3, 
#                 byrow = FALSE, 
#                 species="hg19", 
#                 category = c("chr2"))
add_track(panel_fun = function(gr) {
  # the use of `get_cell_meta_data()` will be introduced later
  chr = get_cell_meta_data("name")  
  grid.rect(gp = gpar(fill = "#EEEEEE"))
  grid.text(chr)
})

add_segments_track(COREtable, COREtable$score, gp = gpar(col = "black", lwd = 4))

fill_string <- NA
if(event == "A"){
  fill_string <- "pink"
} else {
  fill_string <- "blue"
}


add_lines_track(dataInputCORE_density, dataInputCORE_density[[4]], area = TRUE, 
                gp = gpar(fill = fill_string))

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

# TODO: SAVE DATA