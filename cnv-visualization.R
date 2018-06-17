# Download packages
source("https://bioconductor.org/biocLite.R")
biocLite("gtrellis")

# Import libraries
library(gtrellis)
library(circlize)
library(rstudioapi) # load it

cd_local <- function() {
  current_path <- getActiveDocumentContext()$path 
  setwd(dirname(current_path ))
}

cd_doc <- function() {
  setwd("C:/Users/bbece/Documents")
}

cd_local()
samples <- read.table("listLabel.txt", header=T, sep = "\t", stringsAsFactors = F)
for(sample in samples$Organoid){
#sample <- "hN31"
  # Import data
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
  pdf(paste("cnv_image", sample, ".pdf", sep = ""), width = 44, height = 16)
  
  # Visualize bins
  # TODO, switch from segment BED to bin BED
  gtrellis_layout(track_height = c(1,3,3), track_ylim = range(data.frame(facets_bins_data$X.cnlr.)), nrow = 3, n_track = 3, byrow = FALSE, species="hg19")
  add_track(panel_fun = function(gr) {
    # the use of `get_cell_meta_data()` will be introduced later
    chr = get_cell_meta_data("name")  
    grid.rect(gp = gpar(fill = "#EEEEEE"))
    grid.text(chr)
  })
  add_points_track(facets_bins_data, facets_bins_data$X.cnlr., gp = gpar(col = ifelse(facets_bins_data$X.cnlr. > 0.2, "blue", ifelse(facets_bins_data$X.cnlr. > -0.23, "black", "red"))))
  add_segments_track(facets_data, facets_data$X.cnlr.median., gp = gpar(col = ifelse(facets_data$X.cnlr.median. > 0.2, "blue", ifelse(facets_data$X.cnlr.median. > -0.23, "black", "red")), lwd = 4))
  
  # Visualize segments
  #gtrellis_layout(track_ylim= range(data.frame(facets_data$X.cnlr.median.)), nrow = 3, byrow = FALSE)
  # add_segments_track(facets_data, facets_data$X.cnlr.median., gp = gpar(col = ifelse(facets_data$X.cnlr.median. > 0, "red", "green")), lwd = 4)
  #add_segments_track(facets_data, facets_data$X.cnlr.median., gp = gpar(col = ifelse(facets_data$X.cnlr.median. > 0.2, "blue", ifelse(facets_data$X.cnlr.median. > -0.23, "black", "red")), lwd = 4))
  
  dev.off()
}
