args <- commandArgs(trailingOnly = TRUE)

event <- "A" # Should be A (for amplficcation, or D for deletion
outputCsv <- "coreTable.csv"
outputObj <- "newCOREobj.rds"
if (length(args) == 1){
	event <- args[1]
} else if (length(args) == 2){
	event <- args[1]
	outputCsv <- args[2]	
} else if (length(args) == 3){
	event <- args[1]
	outputCsv <- args[2]
	outputObj <- args[3]
}

# Import librariesdd
library(CORE)

print("loaded libraries")

cd_home <- function() {
  setwd("~/code/hN_core_artifacts")
}

cd_home()
samples <- read.table("./resources/sampleList.csv", header=T, sep = "\t", stringsAsFactors = F)
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

print("loaded samples")

chromosomeSizes <- readRDS("./resources/chromosomeSizes.rds")

print("loaded chromosome sizes")

# TODO: Deal with a female XX case (does it matter though?)
rescaleInput <- function(input, chromosomeSizes){
  for(row.index in seq(1, nrow(input))){
    chrom_r <- input[row.index, ]$chrom
    if (is.na(chrom_r) || length(chrom_r) == 0){
      next # TODO: This is to resolve the NA row. Where did it come from?
    } 
    print(chrom_r)
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
  
  loaded_segments <- list(NA)
  loaded_segments.index <- 1
  for(sample in loaded_samples){
    segments <- as.data.frame(read.table(paste("./resources/mappedFacetsFiles/", sample, "--NA12878.mapped.cnv.facets.v0.5.2.bed", sep = ""), header = TRUE, sep=",", stringsAsFactors=FALSE, quote=""))  
    if(event == "A"){
      segments <- segments[segments$X.cnlr. > 0.2,]  
    } else if (event == "D"){
      segments <- segments[segments$X.cnlr. < -0.235,]  
    } else {
      print("No event selected.")
    }
    
    print(head(segments))
    segments <- segments[,c(1, 2, 3)]

    names(segments) <- c("chrom", "start", "end")
    #segments <- rescaleInput(segments, chromosomeSizes)
    
    dataInputCORE <- rbind(dataInputCORE, segments)
  }
  
  # TODO: SKIPPING X AND Y DUE TO INPUT FORMAT ERROR (not accepting string as chr)
  returnme <-  cbind(dataInputCORE)
  returnme$chrom <- as.numeric(returnme$chrom)
  return(returnme)
}

# TODO: Need to verify results - check with Pascal
inputCORE <- generateInputCORE(chromosomeSizes)
print("Printing inputCORE")
print(inputCORE)

print("generated input for CORE")

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

print("generated boundary input for CORE")

#Compute 3 cores and perform no randomization
#(meaningless for estimate of significance).

saveRDS(inputCORE, "inputCORE.rds")

myCOREobj<-CORE(dataIn=inputCORE, maxmark=10, nshuffle=0, seedme=123)

print("Ran first CORE run")
#newCOREobj<-CORE(dataIn=myCOREobj,keep=c("maxmark","seedme","boundaries"),
#                 nshuffle=50,distrib="Rparallel",njobs=4)


newCOREobj<-CORE(dataIn=myCOREobj,keep=c("maxmark","seedme","boundaries"),
 nshuffle=20,distrib="Grid",njobs=2)

print("ran CORE")

saveRDS(newCOREobj, outputObj)

print(paste("Saved CORE obj to", outputObj, sep = ""))

# TODO: Deal with a female XX case (does it matter though?)
rescaleOutput <- function(cores, chromosomeSizes){
  for(row.index in seq(1, nrow(cores))){
    chrom_r <- as.numeric(cores[row.index, ]$chrom)
    print(chrom_r)
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
    #print(total_bp)
    cores[row.index, ]$start <- cores[row.index, ]$start - total_bp
    cores[row.index, ]$end <- cores[row.index, ]$end - total_bp
  }
  return(cores)
}



coreTable <- data.frame(myCOREobj$coreTable)
coreTable <- rescaleOutput(coreTable, chromosomeSizes)
coreTable$chrom <- paste("chr", coreTable$chrom, sep = "")

write.csv(coreTable, outputCsv)

print(paste("Saved coreTable as", outputCsv, sep = ""))


print("done")
