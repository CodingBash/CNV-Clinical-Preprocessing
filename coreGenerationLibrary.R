library(CORE)

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
  
  myCOREobj<- NA
  
  if(missing(inputCOREBoundaries)){
    myCOREobj <- CORE(dataIn=inputCORESegments, maxmark=maxmark, nshuffle=0,
                  seedme=seedme)
  } else {
    myCOREobj <- CORE(dataIn=inputCORESegments, boundaries = inputCOREBoundaries, maxmark=maxmark, nshuffle=0,
                      seedme=seedme)
  }
  
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

selectCnprepSegmentsWithEvent <- function(events, samples, dir, mprob_thresh, bedFormat = TRUE, probes = TRUE, silent = FALSE){
  retrieveCnprepSegmentListFromSamples <- function(dir, samples, silent = FALSE){
    segtables <- data.frame(stringsAsFactors = FALSE)
    for(sample in samples){
      print(paste0(dir, "/", sample, "_segtable.tsv"))
      try({
        segtable.input <- read.table(paste0(dir, "/", sample, "_segtable.tsv"),header=T,as.is=T, stringsAsFactors = FALSE)
        segtables <- rbind(segtables, segtable.input)
      }, silent = silent)
    }
    return(segtables)
  }  
  segmentList <- retrieveCnprepSegmentListFromSamples(dir, samples, silent)
  
  subsetAllCnprepSegmentsByEvent <- function(segmentList, events, mprob_thresh){
    totalSelectedSegments <- data.frame(stringsAsFactors = FALSE)
    
    for(event in events){
      if(event == "A"){
        eventSelectedSegments <- segmentList[segmentList$marginalprob < mprob_thresh & 
                                               segmentList$mediandev > 0,]
        totalSelectedSegments <- rbind(totalSelectedSegments, eventSelectedSegments)
      } else if (event == "D"){
        eventSelectedSegments <- segmentList[segmentList$marginalprob < mprob_thresh & 
                                               segmentList$mediandev < 0,]
        totalSelectedSegments <- rbind(totalSelectedSegments, eventSelectedSegments)
      }
    }
    return(totalSelectedSegments)
  }
  selectedSegments <- subsetAllCnprepSegmentsByEvent(segmentList, events, mprob_thresh)
  
  if(bedFormat == TRUE){
    if(probes == TRUE){
      selectedSegments <- selectedSegments[,c("chrom", "start", "end", "mediandev", "maxzmean", "marginalprob")]
    } else {
      selectedSegments <- selectedSegments[,c("chrom", "abs.pos.start", "abs.pos.end", "mediandev", "maxzmean", "marginalprob")]
    }      
  }
  return(selectedSegments)
}
