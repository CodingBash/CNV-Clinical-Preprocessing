library(CORE)

#
# Given a genome (i.e. hg19), generate the chromosome sizes
#
generateChromosomeSizes <- function(genome){
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

#
# Take a input of multiple segments with the chromosomeSizes, and convert
# segment maploc from chrom.location to absolute.location in bp units
#
chromsomeToAbsoluteBPConversion <- function(input, chromosomeSizes){
  for(row.index in seq(1, nrow(input))){
    chrom_r <- input[row.index, ]$chrom
    if (is.na(chrom_r) || length(chrom_r) == 0){
      next # TODO: This is to resolve the NA row. Where did it come from?
    }
    
    # Calculate BP to add on.
    total_bp <- 0
    if(chrom_r %in% seq(2,22)){ # If chrom is autosomal
      for(i in seq(1, as.numeric(chrom_r) - 1)){ 
        total_bp <- total_bp + chromosomeSizes[paste("chr", i, sep = ""), ]$size # Calculate additional BP from sum of prior autosomal sizes
      }  
    } else if (chrom_r == "X") { # If chrom is X
      for(i in seq(1, 22)){
        total_bp <- total_bp + chromosomeSizes[paste("chr", i, sep = ""), ]$size # Calculate additional BP from sum from all autosomal
      }  
    }  else if (chrom_r == "Y") { # If chrom is Y
      for(i in seq(1, 22)){
        total_bp <- total_bp + chromosomeSizes[paste("chr", i, sep = ""), ]$size # Calculate additional BP from sum from all autosomal and X
      }  
      total_bp <- total_bp + chromosomeSizes["chrX", ]$size
    }
    
    #
    # Update start and end maploc with new absolute location
    #
    input[row.index, ]$start <- input[row.index, ]$start + total_bp
    input[row.index, ]$end <- input[row.index, ]$end + total_bp
  }
  return(input)
}

#
# Take a input of multiple segments (strictly from CORE) with the chromosomeSizes, and convert
# segment maploc from absolute.location to chrom.location in bp units
# 
# TODO: Several duplicate code from vice-versa conversion function. Modularize
# TODO: Allow any bed-file type input instead of just CORE output.
#
absoluteToChromosomeBPConversion <- function(outputCores, chromosomeSizes){
  for(row.index in seq(1, nrow(outputCores))){
    chrom_r <- as.numeric(outputCores[row.index, ]$chrom)
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
    outputCores[row.index, ]$start <- outputCores[row.index, ]$start - total_bp
    outputCores[row.index, ]$end <- outputCores[row.index, ]$end - total_bp
  }
  return(outputCores)
}

#
# Generates the input segments for CORE analysis
#
# @param event - either "A" (amplification) or "D" (deletion) segments extracted
# @param samples - names of samples to retrieve
# @param chromosomeSizes - df of chromosomeSizes for rescaling (if necessary)
# @param dir - directory of the input files
# @param extension - extension of the input files
# TODO: segment filename still too hardcoded. Allow caller to send file name themselves
# @param rescaleInput - indicator if sample segments should be scaled (usually if it is not absolute scaled and is chromsome scaled instead). If TRUE, rescaling with chromosomeSizes is performed
# @param inSampleFolder - ad-hoc solution when the sample file is contained in a folder with sample name
# @param ampCall - lower threshold to call amplifications
# @param delCall - higher threshold to call deletions
#
generateInputCORESegments <- function(event, samples, chromosomeSizes, dir, extension = "cnv.facets.v0.5.2.txt", inSampleFolder = FALSE, rescaleInput = FALSE, ampCall = 0.2, delCall = -0.235){
  dataInputCORE <- data.frame()
  
  loaded_segments <- list(NA)
  loaded_segments.index <- 1
  for(sample in samples){
    segments <- as.data.frame(read.table(paste(dir, if(inSampleFolder == TRUE) paste("Sample_", sample, "/analysis/structural_variants/", sep = "") else "",sample, "--NA12878.", extension, sep = ""), header = TRUE, sep="\t", stringsAsFactors=FALSE, quote=""))
    if(event == "A"){
      segments <- segments[segments$X.cnlr.median. > 0.2,]  
    } else if (event == "D"){
      segments <- segments[segments$X.cnlr.median. < -0.235,]  
    } else {
      print("No event selected.")
    }
    segments <- segments[,c(1, 10, 11)]
    
    names(segments) <- c("chrom", "start", "end")
    
    if(rescaleInput == TRUE){
      segments <- chromsomeToAbsoluteBPConversion(segments, chromosomeSizes)
    }
    
    dataInputCORE <- rbind(dataInputCORE, segments)
  }
  
  returnme <-  cbind(dataInputCORE)
  returnme <- returnme[returnme$chrom != "X" & returnme$chrom != "Y",] # REMOVE X and Y chromosome
  returnme$chrom <- as.numeric(returnme$chrom)
  return(returnme)
}

# 
# Generate BP unit chromosomal boundaries based on chromosomeSizes df
# Similar logic to chromsomeToAbsoluteBPConversion function
#
generateInputCOREBoundaries <- function(chromosomeSizes){
  boundaries <- data.frame(stringsAsFactors = FALSE)
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

#
# Run CORE analysis
#
runCORE <- function(inputCORESegments, inputCOREBoundaries, distrib="Rparallel", maxmark=10, nshuffle=50, seedme, njobs=4) {
  
  myCOREobj<-CORE(dataIn=inputCORESegments, maxmark=maxmark, nshuffle=0,
                  boundaries=inputCOREBoundaries,seedme=seedme)
  
  newCOREobj<-CORE(dataIn=myCOREobj,keep=c("maxmark","seedme","boundaries"),
                   nshuffle=nshuffle,distrib=distrib,njobs=njobs)
  return(newCOREobj)  
}

#
# Format the CORE table from the CORE result object
#
retrieveCORETable <- function(COREobj, chromosomeSizes, rescaleOutput = FALSE){
  COREtable <- data.frame(COREobj$coreTable)
  if(rescaleOutput == TRUE) {
    COREtable <- absoluteToChromosomeBPConversion(COREtable, chromosomeSizes)
  }
  COREtable$chrom <- paste("chr", COREtable$chrom, sep = "")
  return(COREtable)
}
