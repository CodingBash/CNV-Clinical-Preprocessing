#
# Install if needed
#
source("https://bioconductor.org/biocLite.R")
biocLite("Rsamtools")
biocLite("BSgenome.Hsapiens.UCSC.hg19")

#
# Load genome
#
library(BSgenome.Hsapiens.UCSC.hg19)

#
# Load source libraries
# TODO: Organize dependencies
#
source("helperFunctions.R")
source("segmentClusteringLibrary.R")

# TODO: Do not include segments with lower than 5K bp (see paper)

#
# Load input
#
cd_local()
normal_samples <- load_samples(classes = c("N"), sampleList = "sampleList.csv")
cytobands <- retrieveCytobands(dir = "cytoBand.txt")
chromosomeSizes <- generateChromosomeSizes(genome = BSgenome.Hsapiens.UCSC.hg19)
cd_doc()
normalSegments <- selectSegmentsWithEvents(events = c("A", "D", "N"), samples = normal_samples, chromosomeSizes = chromosomeSizes, 
                                              dir = "CSHL/Project_TUV_12995_B01_SOM_Targeted.2018-03-02/", extension = "cnv.facets.v0.5.2.txt", inSampleFolder = TRUE, 
                                              rescaleInput = TRUE, ampCall = 0.2, delCall = -0.235)
tumor_samples <- load_samples(classes = c("T"), sampleList = "sampleList.csv") # TODO: THIS WAS ORIGINAL CLASS "N", RECOMPUTE RESUTS

# Generate norminput argument
norminput <- retrieveNormInput(normalSegments)


for(tumor_samples.i in seq(1, length(tumor_samples))){
  sample <- tumor_samples[tumor_samples.i]
  
  print(paste("Analyzing sample", sample))
  
  #
  # Retrieve sample data
  #
  cd_doc()
  facets_segment_data <- retrieveFacetsSegments(sample, dir = "CSHL/Project_TUV_12995_B01_SOM_Targeted.2018-03-02/")
  facets_snp_data <- retrieveFacetsSnps(sample, dir = "CSHL/Project_TUV_12995_B01_SOM_Targeted.2018-03-02/")
  
  # Generate seginput argument
  seginput <- retrieveSegInput(facets_segment_data, cytobands)
  print(paste("Retrieved segment input for sample", sample))
  
  # Generate ratinput argument
  ratinput <- retrieveRatInput(facets_snp_data)
  print(paste("Retrieved ratio input for sample", sample))
  
  # Run CNprep:CNpreprocessing
  segtable <- runCNpreprocessing(seginput = seginput, ratinput = ratinput, norminput = norminput, modelNames = "E")
  print(paste("Produced segtable for sample", sample))
  
  #
  # Write results out
  #
  cd_local()
  #write.table(segtable, paste("segClusteringResultsPar/", sample, "_segtable.tsv", sep = ""), row.names = F, sep = "\t", quote = FALSE)
  print(head(segtable))
  print(paste("Wrote output for sample", sample))
}

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

chromsomeToAbsoluteBPConversionForSingleEntry <- function(chrom, start, end, chromosomeSizes){
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

#
# Take a input of multiple segments with the chromosomeSizes, and convert
# segment maploc from chrom.location to absolute.location in bp units
#
chromsomeToAbsoluteBPConversion <- function(input, chromosomeSizes){
  for(row.index in seq(1, nrow(input))){
    # Rescale row
    absoluteRow <- chromsomeToAbsoluteBPConversionForSingleEntry(chrom = input[row.index, ]$chrom, start = input[row.index, ]$start, end = input[row.index, ]$end, chromosomeSizes = chromosomeSizes)
    
    #
    # Update start and end maploc with new absolute location
    #
    input[row.index, ]$start <- absoluteRow$start
    input[row.index, ]$end <- absoluteRow$end
  }
  return(input)
}