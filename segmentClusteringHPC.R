library(CNprep)
library(doParallel)

cd_home <- function() {
  setwd("~/code/cnprep_clustering")
}
cd_home()

samples <- read.table("./resources/sampleList.csv", header=T, sep = "\t", stringsAsFactors = F)
cytobands <- read.table("./resources/cytoBand.txt", header=F, sep = "\t", stringsAsFactors = F)
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

for(sample in loaded_samples){
  print(sample)
}

generateInputCORE <- function(chromosomeSizes){
  dataInputCORE <- data.frame()
  
  loaded_segments <- list(NA)
  loaded_segments.index <- 1
  for(sample in loaded_samples){
    segments <- as.data.frame(read.table(paste("./resources/Project_TUV_12995_B01_SOM_Targeted.2018-03-02/Sample_", sample, "/analysis/structural_variants/", sample, "--NA12878.cnv.facets.v0.5.2.txt", sep = ""), header = TRUE, sep="\t", stringsAsFactors=FALSE, quote=""))  
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

chromosomeSizes <- readRDS("./resources/chromosomeSizes.rds")

# TODO: Need to verify results - check with Pascal
normalSegments <- generateInputCORE(chromosomeSizes)



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


classes <- c("T")
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

norminput <- data.frame(stringsAsFactors = FALSE)
for(normalSegments.index in seq(1, nrow(normalSegments))){
  norminput.entry <- data.frame(length = normalSegments[normalSegments.index, ]$X.end. - normalSegments[normalSegments.index, ]$X.start., segmedian = normalSegments[normalSegments.index, ]$X.cnlr.median.)
  norminput <- rbind(norminput, norminput.entry)
}

registerDoParallel()
foreach(loaded_samples.i=1:length(loaded_samples)) %dopar% {
  sample <- loaded_samples[loaded_samples.i]
  
  print(paste("Analyzing sample", sample))
  
  facets_data <- as.data.frame(read.table(paste("./resources/Project_TUV_12995_B01_SOM_Targeted.2018-03-02/Sample_", sample, "/analysis/structural_variants/", sample, "--NA12878.cnv.facets.v0.5.2.txt", sep = ""), header = TRUE, sep="\t", stringsAsFactors=FALSE, quote=""))
  facets_bins_data <- as.data.frame(read.table(paste("./resources/Project_TUV_12995_B01_SOM_Targeted.2018-03-02/Sample_", sample, "/analysis/structural_variants/", sample, "--NA12878.procSample-jseg.cnv.facets.v0.5.2.txt", sep = ""), header = TRUE, sep="\t", stringsAsFactors=FALSE, quote=""))
  
  
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
  
  print(paste("Retrieved segment input for sample", sample))

  ratinput <- data.frame(facets_bins_data$X.cnlr.)
  names(ratinput) <- c(sample)
 
  print(paste("Retrieved ratio input for sample", sample))
  

  
  
  segtable<-CNpreprocessing(segall=seginput,ratall=ratinput,"ID","start","end",
                            chromcol="chrom",bpstartcol="chrom.pos.start",bpendcol="chrom.pos.end",blsize=50,
                            minjoin=0.25,cweight=0.4,bstimes=50,chromrange=1:22,distrib="Rparallel",njobs=40,
                            modelNames="E",normalength=norminput[,1],normalmedian=norminput[,2])
  
  print(paste("Produced segtable for sample", sample))
  
  write.table(segtable, paste("./output/", sample, "_segtable.tsv", sep = ""), row.names = F, sep = "\t", quote = FALSE)
  print(paste("Wrote output for sample", sample))
}
