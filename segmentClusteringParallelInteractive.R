#
# Install if needed
#
#source("https://bioconductor.org/biocLite.R")
#biocLite("Rsamtools")
#biocLite("BSgenome.Hsapiens.UCSC.hg19")

#
# Load genome
#
library(BSgenome.Hsapiens.UCSC.hg19)

#
# Load source libraries
# TODO: Organize dependencies
#
setwd("~/Git-Projects/Git-Research-Projects/drug-response-prediction")
source("helperFunctions.R")
source("segmentClusteringLibrary.R")

# TODO: Do not include segments with lower than 5K bp (see paper)

reference <- "hN30"
res_dir <- "output/FACETS_Reference_hN30_8_2_18_1/"

#
# Load input
#
cd_local()
normal_samples <- load_samples(classes = c("N"), sampleList = "sampleList.csv")
normal_samples <- normal_samples[normal_samples != reference]

cytobands <- retrieveCytobands(dir = "cytoBand.txt")
chromosomeSizes <- generateChromosomeSizes(genome = BSgenome.Hsapiens.UCSC.hg19)

setwd("~/Git-Projects/Git-Research-Projects/FACETS_write_files/")
normalSegments <- selectSegmentsWithEvents(events = c("A", "D", "N"), samples = normal_samples, chromosomeSizes = chromosomeSizes, 
                                           dir = res_dir, sample_subdir="/", reference = reference, extension = "cnv.facets.v0.5.2.txt", inSampleFolder = TRUE, 
                                           rescaleInput = TRUE, ampCall = 0.2, delCall = -0.235)
cd_doc()
# TODO: Does cd_local need to be before this?
tumor_samples <- load_samples(classes = c("T", "F", "M"), sampleList = "sampleList.csv") # TODO: THIS WAS ORIGINAL CLASS "N", RECOMPUTE RESUTS
tumor_samples <- tumor_samples[tumor_samples != reference]
tumor_samples <- tumor_samples[tumor_samples != "hF32"]
tumor_samples <- tumor_samples[tumor_samples != "hF37"]
tumor_samples <- tumor_samples[tumor_samples != "hT105"]

# Generate norminput argument
norminput <- retrieveNormInput(normalSegments)
norminput <- filterNormInput(norminput, length_threshold=10000000) # TODO: Need to be statistically determined

ratinput <- data.frame(stringsAsFactors = FALSE)
seginput <- data.frame(stringsAsFactors = FALSE)
for(tumor_samples.i in seq(1, length(tumor_samples))){
  sample <- tumor_samples[tumor_samples.i]
  
  print(paste("Analyzing sample", sample))
  
  #
  # Retrieve sample data
  #
  setwd("~/Git-Projects/Git-Research-Projects/FACETS_write_files/")
  facets_segment_data <- retrieveFacetsSegments(sample, sample_subdir = "/", reference = reference, dir = res_dir)
  facets_snp_data <- retrieveFacetsSnps(sample, sample_subdir = "/", reference = reference, dir = res_dir)
  
  # Generate seginput argument
  seginput.entry <- retrieveSegInput(facets_segment_data, sample, chromosomeSizes, cytobands)
  print(paste("Retrieved segment input for sample", sample))

  seginput <- rbind(seginput, seginput.entry)
  
  # Generate ratinput argument
  ratinput.entry <- retrieveRatInput(facets_snp_data, sample)
  print(paste("Retrieved ratio input for sample", sample))
  
  if(dim(ratinput)[1] == 0){
    ratinput <- cbind(ratinput.entry)
  } else{
    ratinput <- cbind(ratinput, ratinput.entry)
  }
}

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