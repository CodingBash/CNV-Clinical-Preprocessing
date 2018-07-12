# Download packages
source("https://bioconductor.org/biocLite.R")
biocLite("gtrellis")
biocLite("ComplexHeatmap")

# Import libraries
library(gtrellis)
library(circlize)
library(rstudioapi) # load it
library(ComplexHeatmap)

cd_local <- function() {
  current_path <- getActiveDocumentContext()$path 
  setwd(dirname(current_path ))
}

cd_doc <- function() {
  setwd("C:/Users/bbece/Documents")
}

cd_local()
samples <- read.table("sampleList.csv", header=T, sep = "\t", stringsAsFactors = F)
#for(sample in samples$Organoids){
sample <- "hT30"
# Import data
segment_clusters <- as.data.frame(read.table(paste("segClusteringResults/", sample, "_segtable.tsv", sep = ""), header = TRUE, sep="\t", stringsAsFactors=FALSE, quote=""))
segment_clusters <- segment_clusters[,c(6,7,8,20)]
segment_clusters <- segment_clusters[segment_clusters$chrom != "X",]
segment_clusters$chrom <- paste("chr", segment_clusters$chrom, sep = "")


cd_doc()
#facets_data <- as.data.frame(read.table("CSHL/Project_TUV_12995_B01_SOM_Targeted.2018-03-02/Sample_hT1/analysis/structural_variants/hT1--NA12878.cnv.facets.v0.5.2.txt", header = TRUE, sep="\t", stringsAsFactors=FALSE, quote=""))
facets_data <- as.data.frame(read.table(paste("CSHL/Project_TUV_12995_B01_SOM_Targeted.2018-03-02/Sample_", sample, "/analysis/structural_variants/", sample, "--NA12878.cnv.facets.v0.5.2.txt", sep = ""), header = TRUE, sep="\t", stringsAsFactors=FALSE, quote=""))
facets_bins_data <- as.data.frame(read.table(paste("CSHL/Project_TUV_12995_B01_SOM_Targeted.2018-03-02/Sample_", sample, "/analysis/structural_variants/", sample, "--NA12878.procSample-jseg.cnv.facets.v0.5.2.txt", sep = ""), header = TRUE, sep="\t", stringsAsFactors=FALSE, quote=""))
# data_bed <- as.data.frame(read.table(paste("CSHL/Project_TUV_12995_B01_SOM_Targeted.2018-03-02/additional_analyses/cnv/h", tissue_type, tissue_number, "--NA12878.facets.bed", sep = ""), header = FALSE, sep="\t", stringsAsFactors=FALSE, quote=""))
#data_bed$V1 <- paste("chr", data_bed$V1, sep="")


facets_data <- facets_data[,c(1, 10, 11, 5)]
facets_data$X.chrom. <- paste("chr", facets_data$X.chrom., sep="")

facets_bins_data$start <- seq(1, length.out=nrow(facets_bins_data), by=1)
facets_bins_data$end <- seq(2, length.out=nrow(facets_bins_data), by=1)
facets_bins_data <- facets_bins_data[,c(1, 2, 2, 11)]
facets_bins_data$X.chrom. <- paste("chr", facets_bins_data$X.chrom., sep="")


print(facets_data)
print(facets_bins_data)
# print(data_bed)


cd_local()
#pdf(paste("cnv_image_merged_chr19_", sample, ".pdf", sep = ""), width = 44, height = 16)

# Visualize bins
# TODO, switch from segment BED to bin BED
lgd = Legend(at = c("duplication", "nuetral", "deletion"), title = "Class", type = "lines", legend_gp = gpar(col = c("orange", "blue", "red")))
gtrellis_layout(track_height = c(1,8,1),
                track_axis = c(FALSE, TRUE, FALSE), 
                track_ylim = range(seq(-6,2)), 
                nrow = 3, 
                n_track = 3, 
                byrow = FALSE, 
                species="hg19",
                # category = "chr8",
                legend = lgd)

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
add_points_track(facets_bins_data, facets_bins_data$X.cnlr., gp = gpar(fill = "gray"))
#add_points_track(facets_bins_data, facets_bins_data$X.cnlr., gp = gpar(col = ifelse(facets_bins_data$X.cnlr. > 0.2, "blue", ifelse(facets_bins_data$X.cnlr. > -0.23, "black", "red"))))
add_segments_track(segment_clusters, segment_clusters$maxzmean, track = current_track(), gp = gpar(col = ifelse(segment_clusters$maxzmean > 0.2, "orange", ifelse(segment_clusters$maxzmean > -0.23, "blue", "red")), lwd = 4))

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
# Visualize segments
#gtrellis_layout(track_ylim= range(data.frame(facets_data$X.cnlr.median.)), nrow = 3, byrow = FALSE)
# add_segments_track(facets_data, facets_data$X.cnlr.median., gp = gpar(col = ifelse(facets_data$X.cnlr.median. > 0, "red", "green")), lwd = 4)
#add_segments_track(facets_data, facets_data$X.cnlr.median., gp = gpar(col = ifelse(facets_data$X.cnlr.median. > 0.2, "blue", ifelse(facets_data$X.cnlr.median. > -0.23, "black", "red")), lwd = 4))

#dev.off()
#}
