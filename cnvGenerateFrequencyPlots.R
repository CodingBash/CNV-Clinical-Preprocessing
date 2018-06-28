# Download packages
source("https://bioconductor.org/biocLite.R")
biocLite("gtrellis")
biocLite("ComplexHeatmap")

# Import libraries
library(gtrellis)
library(circlize)
library(rstudioapi) # load it
library(ComplexHeatmap)
library(RColorBrewer)

cd_local <- function() {
  current_path <- getActiveDocumentContext()$path 
  setwd(dirname(current_path ))
}

cd_doc <- function() {
  setwd("C:/Users/bbece/Documents")
}

cd_local()
samples <- read.table("sampleList.csv", header=T, sep = "\t", stringsAsFactors = F)
classes <- c("N")
event <- "A" # A - amplification, D - deletion

loaded_samples <- c(NA)
loaded_samples.index <- 1
for(sample in samples$Organoids){
  for(class in classes){
    if(substring(sample, 2,2) == class){
      loaded_samples[loaded_samples.index] <- sample
      loaded_samples.index <- loaded_samples.index + 1
      next
    }
  }
}


  ########### TODO: Make into method (also this piece of code is modelled after cnvCores.R#generateInputCORE())
  dataInputCORE <- data.frame()
  
  cd_doc()
  loaded_bins <- list(NA)
  loaded_bins.index <- 1
  for(sample in loaded_samples){
    bins <- as.data.frame(read.table(paste("CSHL/Project_TUV_12995_B01_SOM_Targeted.2018-03-02/Sample_", sample, "/analysis/structural_variants/", sample, "--NA12878.procSample-jseg.cnv.facets.v0.5.2.txt", sep = ""), header = TRUE, sep="\t", stringsAsFactors=FALSE, quote=""))
    head(bins)
    if(event == "A"){
      bins <- bins[bins$X.cnlr. > 0.2,]  
    } else if (event == "D"){
      bins <- bins[bins$X.cnlr. < -0.235,]  
    }
    
    bins <- bins[,c(1, 2, 2, 11)]
    names(bins) <- c("chrom", "start", "end", "cnlr")
    #segments <- rescaleInput(segments, chromosomeSizes)
    
    dataInputCORE <- rbind(dataInputCORE, bins)
  }
  
  # TODO: SKIPPING X AND Y DUE TO INPUT FORMAT ERROR (not accepting string as chr)
  dataInputCORE <- dataInputCORE[dataInputCORE$chrom != "X" & dataInputCORE$chrom != "Y",]
  dataInputCORE$chrom <- paste("chr", dataInputCORE$chrom, sep = "")
  
  ##########
  dataInputCORE_density = circlize::genomicDensity(dataInputCORE, window.size = 5e6)
  head(dataInputCORE_density)
  
  
  
  
  gtrellis_layout(n_track = 4, nrow = 3, byrow = FALSE,
                  track_axis = c(FALSE, TRUE, FALSE, FALSE), 
                  track_height = unit.c(2*grobHeight(textGrob("chr1")), 
                                        unit(1.5, "null"), 
                                        unit(5, "mm"), 
                                        unit(3, "mm")), 
                  track_ylim = c(0, 1, c(0, max(dataInputCORE_density[[4]])), 0, 1, 0, 1),
                  track_ylab = c("", "density", "", ""))
  
  # track for chromosome names
  add_track(panel_fun = function(gr) {
    # the use of `get_cell_meta_data()` will be introduced later
    chr = get_cell_meta_data("name")  
    grid.rect(gp = gpar(fill = "#EEEEEE"))
    grid.text(chr)
  })
  
  fill_string <- NA
  if(event == "A"){
    fill_string <- "pink"
  } else {
    fill_string <- "blue"
  }
  
  
  # track for genomic density
  add_lines_track(dataInputCORE_density, dataInputCORE_density[[4]], area = TRUE, 
                  gp = gpar(fill = fill_string))
  
  # TODO: Set default to lowest heat
  col_fun = circlize::colorRamp2(seq(0, max(dataInputCORE_density[[4]]), length = 11), 
                                 rev(brewer.pal(11, "RdYlBu")))
  add_heatmap_track(dataInputCORE_density, dataInputCORE_density[[4]], fill = col_fun)
  
  # track for ideogram
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

