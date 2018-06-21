source("https://bioconductor.org/biocLite.R")
biocLite("BSgenome.Hsapiens.UCSC.hg19")

install.packages("CORE")
install.packages("rstudioapi")
# Import libraries
library(rstudioapi) # load it
library(BSgenome.Hsapiens.UCSC.hg19)
library(CORE)
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

# A table of chromosome boundary positions for DNA copy number analysis
generateChromosomeSizes <- function(){
  genome <- BSgenome.Hsapiens.UCSC.hg19
  seqlengths(genome) <- seqlengths(Hsapiens)
  
  # Create chromosome vector
  chrom_vec <- c(NA)
  chrom_vec.index <- 1
  for(i in append(seq(1,22, by=1), c("X", "Y"))){
    chrom_vec[chrom_vec.index] <- paste("chr", i, sep = "")  
    chrom_vec.index <- chrom_vec.index + 1
  }
  chromosomeSizes <- data.frame()
  
  for(chrom_i in chrom_vec){
    df = data.frame(chrom = chrom_i, size = seqlengths(genome)[chrom_i])
    chromosomeSizes <- rbind(chromosomeSizes, df)
  }
  return(chromosomeSizes)
}

chromosomeSizes <- generateChromosomeSizes()


# TODO: Deal with a female XX case (does it matter though?)
rescaleInput <- function(input, chromosomeSizes){
  for(row.index in seq(1, nrow(input))){
    chrom_r <- input[row.index, ]$chrom
    total_bp <- 0
    if(chrom_r %in% seq(2,22)){
      for(i in seq(1, as.numeric(chrom_r) - 1)){
        total_bp <- total_bp + chromosomeSizes[paste("chr", i, sep = ""), ]$size
      }  
    } else if (chrom_r == "X") {
      for(i in seq(1, 22)){
        total_bp <- total_bp + chromosomeSizes[paste("chr", i, sep = ""), ]$size
      }  
    }  else if (chrom_r == "Y") {
      for(i in seq(1, 22)){
        total_bp <- total_bp + chromosomeSizes[paste("chr", i, sep = ""), ]$size
      }  
      total_bp <- total_bp + chromosomeSizes["chrX", ]$size
    }
    input[row.index, ]$start <- input[row.index, ]$start + total_bp
    input[row.index, ]$end <- input[row.index, ]$end + total_bp
  }
  return(input)
}

# A table of DNA copy number gain events observed in 100 individual tumor cells
generateInputCORE <- function(chromosomeSizes){
  dataInputCORE <- data.frame()
  
  cd_doc()
  loaded_segments <- list(NA)
  loaded_segments.index <- 1
  for(sample in loaded_samples){
    segments <- as.data.frame(read.table(paste("CSHL/Project_TUV_12995_B01_SOM_Targeted.2018-03-02/Sample_", sample, "/analysis/structural_variants/", sample, "--NA12878.cnv.facets.v0.5.2.txt", sep = ""), header = TRUE, sep="\t", stringsAsFactors=FALSE, quote=""))  
    segments <- segments[,c(1, 10, 11)]
    names(segments) <- c("chrom", "start", "end")
    segments <- rescaleInput(segments, chromosomeSizes)
    
    dataInputCORE <- rbind(dataInputCORE, segments)
  }
  
  # TODO: SKIPPING X AND Y DUE TO INPUT FORMAT ERROR (not accepting string as chr)
  returnme <- dataInputCORE[dataInputCORE$chrom != "X" & dataInputCORE$chrom != "Y",]
  returnme$chrom <- as.numeric(returnme$chrom)
  return(returnme)
}

# TODO: Need to verify results - check with Pascal
inputCORE <- generateInputCORE(chromosomeSizes)

rescaleBoundaries <- function(chromosomeSizes){
  boundaries <- data.frame(stringsAsFactors = FALSE)
  # TODO: This may not work, may need to paste "chr" --- it works fine
  for(row.index in seq(1, nrow(chromosomeSizes))){
    chrom_r <- chromosomeSizes[row.index, ]$chrom
    total_bp <- 0
    last_size <- chromosomeSizes[chrom_r, ]$size
    
    if(chrom_r %in% paste("chr", seq(1,22), sep = "")){ # TODO: chrom_r index may be off
      for(i in paste("chr", seq(1, as.numeric(substring(chrom_r, 4))), sep = "")){
        total_bp <- total_bp + chromosomeSizes[i, ]$size
      }  
    } else if (chrom_r == "chrX") {
      for(i in paste("chr", seq(1, 22), sep = "")){
        total_bp <- total_bp + chromosomeSizes[i, ]$size
      }  
    }  else if (chrom_r == "chrY") {
      for(i in paste("chr", seq(1, 22), sep = "")){
        total_bp <- total_bp + chromosomeSizes[i, ]$size
      }  
      total_bp <- total_bp + chromosomeSizes["chrX", ]$size
    } else {
      next
    }
    df = data.frame(chrom = substring(chrom_r, 4), start = total_bp - last_size + 1, end = total_bp, stringsAsFactors = FALSE)
    boundaries <- rbind(boundaries, df)
  }
  
  # TODO: SKIPPING X AND Y DUE TO INPUT FORMAT ERROR (not accepting string as chr)
  returnme <- boundaries[boundaries$chrom != "X" & boundaries$chrom != "Y",]
  returnme$chrom <- as.numeric(returnme$chrom)
  
  return(returnme)
}
inputBoundaries <- rescaleBoundaries(chromosomeSizes)


#Compute 3 cores and perform no randomization
#(meaningless for estimate of significance).

myCOREobj<-CORE(dataIn=inputCORE, maxmark=10, nshuffle=0,
                boundaries=inputBoundaries,seedme=123)
## Not run:
#Extend this computation to a much larger number of randomizations,
#using 2 cores of a host computer.
newCOREobj<-CORE(dataIn=myCOREobj,keep=c("maxmark","seedme","boundaries"),
                 nshuffle=50,distrib="Rparallel",njobs=4)
#When using "Grid", make sure you have write
## Not run:
#Extend this computation to a much larger number of randomizations,
#using 2 cores of a host computer.
# newCOREobj<-CORE(dataIn=myCOREobj,keep=c("maxmark","seedme","boundaries"),
# nshuffle=20,distrib="Rparallel",njobs=2)
#When using "Grid", make sure you have write premission to the current
#work space.
# newCOREobj<-CORE(dataIn=myCOREobj,keep=c("maxmark","seedme","boundaries"),
# nshuffle=20,distrib="Grid",njobs=2)
## End(Not run)

install.packages("devtools")
library(devtools)
install_github("wefang/ghelper")

boundaries <- data.frame(myCOREobj$coreTable)
boundaries$chrom <- paste("chr", boundaries$chrom, sep = "")

library(gtrellis)
library(circlize)
library(rstudioapi) # load it
library(ComplexHeatmap)

lgd = Legend(at = c("duplication", "nuetral", "deletion"), title = "Class", type = "lines", legend_gp = gpar(col = c("orange", "blue", "red")))
gtrellis_layout(track_height = c(2,5,1),
                track_axis = c(FALSE, TRUE, FALSE), 
                track_ylim = range(seq(-1,1)), 
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

add_segments_track(boundaries, 0, gp = gpar(col = "black", lwd = 4))

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