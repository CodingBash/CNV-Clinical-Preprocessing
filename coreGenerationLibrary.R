library(CORE)

#
# Take a input of multiple segments (strictly from CORE) with the chromosomeSizes, and convert
# segment maploc from absolute.location to chrom.location in bp units
# 
# TODO: Several duplicate code from vice-versa conversion function. Modularize
# TODO: Allow any bed-file type input instead of just CORE output, then move to helperFunction.R script
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
