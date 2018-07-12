install.packages("rstudioapi")

# Import libraries
library(rstudioapi) # load it
# Import libraries
library(gtrellis)
library(circlize)
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
classes <- c("N")

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


# A table of DNA copy number gain events observed in 100 individual tumor cells
generateInputCORE <- function(){
  dataInputCORE <- data.frame()
  
  cd_doc()
  loaded_segments <- list(NA)
  loaded_segments.index <- 1
  for(sample in loaded_samples){
    segments <- as.data.frame(read.table(paste("CSHL/Project_TUV_12995_B01_SOM_Targeted.2018-03-02/Sample_", sample, "/analysis/structural_variants/", sample, "--NA12878.cnv.facets.v0.5.2.txt", sep = ""), header = TRUE, sep="\t", stringsAsFactors=FALSE, quote=""))  
    #print(segments)
    segments <- segments[,c(1, 10, 11, 5)]
    #print(segments)
    names(segments) <- c("chrom", "start", "end", "cnlr")
    dataInputCORE <- rbind(dataInputCORE, segments)
  }
  dataInputCORE$chrom <- paste("chr", dataInputCORE$chrom, sep = "")
  return(dataInputCORE)
}

testInputCORE <- generateInputCORE()

lgd = Legend(at = c("duplication", "nuetral", "deletion"), title = "Class", type = "lines", legend_gp = gpar(col = c("orange", "blue", "red")))
gtrellis_layout(track_height = c(2,5,1),
                track_axis = c(FALSE, TRUE, FALSE), 
                track_ylim = range(data.frame(testInputCORE$cnlr)), 
                nrow = 3, 
                n_track = 3, 
                byrow = FALSE, 
                species="hg19",
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
#add_points_track(facets_bins_data, facets_bins_data$X.cnlr., gp = gpar(col = ifelse(facets_bins_data$X.cnlr. > 0.2, "blue", ifelse(facets_bins_data$X.cnlr. > -0.23, "black", "red"))))
add_segments_track(testInputCORE, testInputCORE$cnlr, gp = gpar(col = ifelse(testInputCORE$cnlr > 0.2, "orange", ifelse(testInputCORE$cnlr > -0.23, "blue", "red")), lwd = 4))

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

