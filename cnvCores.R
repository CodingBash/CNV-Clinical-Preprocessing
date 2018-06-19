source("https://bioconductor.org/biocLite.R")
biocLite("BSgenome.Hsapiens.UCSC.hg19")

install.packages("CORE")
install.packages("rstudioapi")
# Import libraries
library(rstudioapi) # load it
library(BSgenome.Hsapiens.UCSC.hg19)

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
    segments <- segments[,c(1, 10, 11)]
    names(segments) <- c("chrom", "start", "end")
    dataInputCORE <- rbind(dataInputCORE, segments)
  }
  return(dataInputCORE)
}

# A table of chromosome boundary positions for DNA copy number analysis
generateInputBoundaries <- function(){
  genome <- BSgenome.Hsapiens.UCSC.hg19
  seqlengths(genome) <- seqlengths(Hsapiens)
  
  # Create chromosome vector
  chrom_vec <- c(NA)
  chrom_vec.index <- 1
  for(i in append(seq(1,22, by=1), c("X", "Y"))){
    chrom_vec[chrom_vec.index] <- paste("chr", i, sep = "")  
    chrom_vec.index <- chrom_vec.index + 1
  }
  dataInputBoundaries <- data.frame()
  
  for(chrom_i in chrom_vec){
    df = data.frame(chrom = chrom_i, start = 0, end = seqlengths(genome)[chrom_i])
    dataInputBoundaries <- rbind(dataInputBoundaries, df)
  }
  return(dataInputBoundaries)
}

testInputCORE <- generateInputCORE()
testInputBoundaries <- generateInputBoundaries()

#Compute 3 cores and perform no randomization
#(meaningless for estimate of significance).
data(testInputCORE)
data(testInputBoundaries)
myCOREobj<-CORE(dataIn=testInputCORE,maxmark=3,nshuffle=0,
boundaries=testInputBoundaries,seedme=123)
## Not run:
#Extend this computation to a much larger number of randomizations,
#using 2 cores of a host computer.
newCOREobj<-CORE(dataIn=myCOREobj,keep=c("maxmark","seedme","boundaries"),
nshuffle=20,distrib="Rparallel",njobs=2)
#When using "Grid", make sure you have write premission to the current
#work space.
newCOREobj<-CORE(dataIn=myCOREobj,keep=c("maxmark","seedme","boundaries"),
nshuffle=20,distrib="Grid",njobs=2)
## End(Not run)
