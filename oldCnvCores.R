#
# This script generates CORES from input files in BP units (instead of SNP/bin units).
# This script is interactive and is NOT meant to be ran on a HPC job from a shell script.
#

#
# Only run if not installed
#
source("https://bioconductor.org/biocLite.R")
biocLite("BSgenome.Hsapiens.UCSC.hg19")
install.packages("CORE")

library(BSgenome.Hsapiens.UCSC.hg19)
library(CORE)

source("coreGenerationLibrary.R")
source("helperFunctions.R")

event <- "D" # "A" - amplification, "D" - deletion

cd_local()
samples <- load_samples(classes = c("T"), sampleList = "sampleList.csv")
chromosomeSizes <- generateChromosomeSizes(genome = BSgenome.Hsapiens.UCSC.hg19)
cd_doc()
inputCORESegments <- generateInputCORESegments(event, samples, chromosomeSizes, dir = "CSHL/Project_TUV_12995_B01_SOM_Targeted.2018-03-02/", extension = "cnv.facets.v0.5.2.txt", inSampleFolder = TRUE, rescaleInput = TRUE, ampCall = 0.2, delCall = -0.235)
inputCOREBoundaries <- generateInputCOREBoundaries(chromosomeSizes)

outputCOREobj <- runCORE(inputCORE, inputCOREBoundaries, distrib="Rparallel", maxmark=10, nshuffle=20, seedme=123, njobs=4)

retrieveCORETable(outputCOREobj, chromosomeSizes, rescaleOutput = TRUE)

# TODO: Verify outputCOREobj results, and create Jupyter notebook for CORE visualization

lgd = Legend(at = c("duplication", "nuetral", "deletion"), title = "Class", type = "lines", legend_gp = gpar(col = c("orange", "blue", "red")))
gtrellis_layout(track_height = unit.c(2*grobHeight(textGrob("chr1")), 
                                      unit(1.0, "null"), 
                                      unit(0.5, "null"), 
                                      unit(3, "mm")),
                track_axis = c(FALSE, TRUE, TRUE, FALSE), 
                track_ylim = c(0, 1, range(data.frame(coreTable$score)), c(0, max(dataInputCORE_density[[4]])), 0, 1),
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

add_segments_track(coreTable, coreTable$score, gp = gpar(col = "black", lwd = 4))

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

