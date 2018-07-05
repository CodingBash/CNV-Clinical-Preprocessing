source("https://bioconductor.org/biocLite.R")
biocLite("BSgenome.Hsapiens.UCSC.hg19")

install.packages("CNprep")
library(BSgenome.Hsapiens.UCSC.hg19)
library(CNprep)
library(rstudioapi) # load it


# TODO: Remove the X and Y chromosome segments/bins
# TODO: Do not include segments with lower than 5K bp (see paper)

# TODO: A lot of the code below is duplicate between many of my R scripts. Need to modularize
cd_local <- function() {
  current_path <- getActiveDocumentContext()$path 
  setwd(dirname(current_path ))
}

cd_doc <- function() {
  setwd("C:/Users/bbece/Documents")
}

cd_local()
samples <- read.table("sampleList.csv", header=T, sep = "\t", stringsAsFactors = F)
cytobands <- read.table("cytoBand.txt", header=F, sep = "\t", stringsAsFactors = F)
names(cytobands) <- c("chrom", "start", "end", "cytoloc", "stain")

# TODO: This is a duplicate from cnvCores.R!!! Modularize
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

generateInputCORE <- function(chromosomeSizes){
  dataInputCORE <- data.frame()
  
  cd_doc()
  loaded_segments <- list(NA)
  loaded_segments.index <- 1
  for(sample in loaded_samples){
    segments <- as.data.frame(read.table(paste("CSHL/Project_TUV_12995_B01_SOM_Targeted.2018-03-02/Sample_", sample, "/analysis/structural_variants/", sample, "--NA12878.cnv.facets.v0.5.2.txt", sep = ""), header = TRUE, sep="\t", stringsAsFactors=FALSE, quote=""))  
    # TODO: CHANGED from cnvCores.R - removed event filtering
    #if(event == "A"){
    #  segments <- segments[segments$X.cnlr.median. > 0.2,]  
    #} else if (event == "D"){
    #  segments <- segments[segments$X.cnlr.median. < -0.235,]  
    #}
    
    #segments <- segments[,c(1, 10, 11)] TODO: CHANGED FROM cnvCores.R version
    #names(segments) <- c("chrom", "start", "end") TODO: CHANGED FROM cnvCores.R version
    #segments <- rescaleInput(segments, chromosomeSizes)
    
    dataInputCORE <- rbind(dataInputCORE, segments)
  }
  
  # TODO: SKIPPING X AND Y DUE TO INPUT FORMAT ERROR (not accepting string as chr)                          IMPORTANT
  returnme <- dataInputCORE[dataInputCORE$X.chrom. != "X" & dataInputCORE$X.chrom. != "Y",] # TODO: CHANGED FROM cnvCores.R -> $chrom to #X.chrom.
  
  # TODO: CHANGED FROM cnvCores.R version ... intentionally skipping X and Y chromosome
  #returnme <- cbind(dataInputCORE)
  #if(length(returnme[returnme$chrom == 'X',]$chrom) > 0) returnme[returnme$chrom == 'X',]$chrom <- "23"
  #if(length(returnme[returnme$chrom == 'Y',]$chrom) > 0) returnme[returnme$chrom == 'Y',]$chrom <- "24"
  returnme$X.chrom. <- as.numeric(returnme$X.chrom.) # TODO: CHANGED from cnvCores.R -> changed from $chrom to $X.chrom.
  return(returnme)
}

# TODO: Need to verify results - check with Pascal
normalSegments <- generateInputCORE(chromosomeSizes)



cd_doc()
facets_data <- as.data.frame(read.table(paste("CSHL/Project_TUV_12995_B01_SOM_Targeted.2018-03-02/Sample_", sample, "/analysis/structural_variants/", sample, "--NA12878.cnv.facets.v0.5.2.txt", sep = ""), header = TRUE, sep="\t", stringsAsFactors=FALSE, quote=""))
facets_bins_data <- as.data.frame(read.table(paste("CSHL/Project_TUV_12995_B01_SOM_Targeted.2018-03-02/Sample_", sample, "/analysis/structural_variants/", sample, "--NA12878.procSample-jseg.cnv.facets.v0.5.2.txt", sep = ""), header = TRUE, sep="\t", stringsAsFactors=FALSE, quote=""))


#facets_data <- facets_data[,c(1, 10, 11, 5)]
#facets_data$X.chrom. <- paste("chr", facets_data$X.chrom., sep="")

#facets_bins_data$start <- seq(1, length.out=nrow(facets_bins_data), by=1)
#facets_bins_data$end <- seq(2, length.out=nrow(facets_bins_data), by=1)
#facets_bins_data <- facets_bins_data[,c(1, 2, 2, 11)]
#facets_bins_data$X.chrom. <- paste("chr", facets_bins_data$X.chrom., sep="")


head(facets_data)
head(facets_bins_data)
# print(data_bed)

# TODO: Duplicate CODE from cnvCores.R! Need to modularize
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


# TODO: Duplicate code! -> The only thing changed is pass in a entry's chrom, start, and end instead of whole input
# --- so we remove the input iteratation. Also, the result is stored in a 1-row dataframe and returned. Essentially, this is the same method but
# --- but works for just a single row. We can modularize by having the method that works on a whole input just utilize this one.
# TODO: Deal with a female XX case (does it matter though?)
rescaleInput <- function(chrom, start, end, chromosomeSizes){
    chrom_r <- chrom
    if (is.na(chrom_r) || length(chrom_r) == 0){
      next # TODO: This is to resolve the NA row. Where did it come from?
    } 
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
    abs_start <- start + total_bp
    abs_end <- end + total_bp
    returnme <- data.frame(start = abs_start, end = abs_end)
    return(returnme)
}

findCytolocation <- function(chrom, chrom.position){
  row <- cytobands[cytobands$chrom == paste("chr", chrom, sep = "") & cytobands$start <= chrom.position & cytobands$end >= chrom.position, ]
  returnme <- data.frame(row)$cytoloc
  return(returnme)
}

sample <- "hT30"
# GENERATE SEGINPUT FROM FACETS SAMPLE
seginput <- data.frame(stringsAsFactors = FALSE)
for(facets_data.index in seq(1, nrow(facets_data))){
  abs_position <- rescaleInput(facets_data[facets_data.index,]$X.chrom., facets_data[facets_data.index,]$X.start., facets_data[facets_data.index,]$X.end., chromosomeSizes)
  
  probes.start = 0
  probes.end = 0
  for(i in seq(1, facets_data.index)){
    if(i != facets_data.index){
      probes.start <- probes.start + facets_data[i,]$X.num.mark.
    }
    probes.end <- probes.end + facets_data[i,]$X.num.mark.
  }
  probes.start <- probes.start + 1
  
  
  cytoband.my.start <- findCytolocation(chrom = facets_data[facets_data.index,]$X.chrom., chrom.position = facets_data[facets_data.index,]$X.start.)
  cytoband.my.end <- findCytolocation(chrom = facets_data[facets_data.index,]$X.chrom., chrom.position = facets_data[facets_data.index,]$X.end.)
  
  
  seginput.entry <- data.frame(ID = sample, start = probes.start, end = probes.end, 
                               num.probes = facets_data[facets_data.index,]$X.num.mark., seg.median = facets_data[facets_data.index,]$X.cnlr.median., 
                               chrom = facets_data[facets_data.index,]$X.chrom., chrom.pos.start = facets_data[facets_data.index,]$X.start., 
                               chrom.pos.end = facets_data[facets_data.index,]$X.end., cytoband.start = cytoband.my.start, 
                               cytoband.end = cytoband.my.end, abs.pos.start = abs_position$start,
                               abs.pos.end = abs_position$end)
  
  
  seginput <- rbind(seginput, seginput.entry)
}

ratinput <- data.frame(stringsAsFactors = FALSE)
for(facets_bins_data.index in seq(1, nrow(facets_bins_data))){
  ratinput.entry <- data.frame(sample = facets_bins_data[facets_bins_data.index,]$X.cnlr.)
  names(ratinput.entry) <- c(sample)
  ratinput <- rbind(ratinput, ratinput.entry)
}

norminput <- data.frame(stringsAsFactors = FALSE)
for(normalSegments.index in seq(1, nrow(normalSegments))){
  norminput.entry <- data.frame(length = normalSegments[normalSegments.index, ]$X.end. - normalSegments[normalSegments.index, ]$X.start., segmedian = normalSegments[normalSegments.index, ]$X.cnlr.median.)
  norminput <- rbind(norminput, norminput.entry)
}



segtable<-CNpreprocessing(segall=seginput,ratall=ratinput,"ID","start","end",
                          chromcol="chrom",bpstartcol="chrom.pos.start",bpendcol="chrom.pos.end",blsize=50,
                          minjoin=0.25,cweight=0.4,bstimes=50,chromrange=1:22,distrib="Rparallel",njobs=40,
                          modelNames="E",normalength=norminput[,1],normalmedian=norminput[,2])


###############
data(segexample)
data(ratexample)
data(normsegs)


#small toy example
segtable<-CNpreprocessing(segall=segexample[segexample[,"ID"]=="WZ1",],
                          ratall=ratexample,"ID","start","end",chromcol="chrom",bpstartcol="chrom.pos.start",
                          bpendcol="chrom.pos.end",blsize=50,minjoin=0.25,cweight=0.4,bstimes=50,
                          chromrange=1:3,distrib="Rparallel",njobs=2,modelNames="E",
                          normalength=normsegs[,1],normalmedian=normsegs[,2])
## Not run:
#Example 1: 5 whole genome analysis, choosing the right format of arguments
segtable<-CNpreprocessing(segall=segexample,ratall=ratexample,"ID","start","end",
                          chromcol="chrom",bpstartcol="chrom.pos.start",bpendcol="chrom.pos.end",blsize=50,
                          minjoin=0.25,cweight=0.4,bstimes=50,chromrange=1:22,distrib="Rparallel",njobs=40,
                          modelNames="E",normalength=normsegs[,1],normalmedian=normsegs[,2])
#Example 2: how to use annotexample, when segment table does not have columns of
#integer postions in terms of measuring units(probes), such as "mysegs" below
mysegs<-segexample[,c(1,5:12)]
data(annotexample)
segtable<-CNpreprocessing(segall=mysegs,ratall=ratexample,"ID",chromcol="chrom",
                          bpstartcol="chrom.pos.start",bpendcol="chrom.pos.end",annot=annotexample,
                          annotstartcol="CHROM.POS",annotendcol="CHROM.POS",annotchromcol="CHROM",
                          blsize=50,minjoin=0.25,cweight=0.4,bstimes=50,chromrange=1:22,distrib="Rparallel",
                          njobs=40,modelNames="E",normalength=normsegs[,1],normalmedian=normsegs[,2])
## End(Not run)

data(ratexample)
#Plot the whole genome log ratio data for the first profile
#Note X and Y chromosomes at the far right of the plot
plot(ratexample[,1])
