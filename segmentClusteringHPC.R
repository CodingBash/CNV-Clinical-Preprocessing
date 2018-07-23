#
# Load source libraries
# TODO: Organize dependencies
#
setwd("~/code/cnprep_clustering/scripts")
source("helperFunctions.R")
source("segmentClusteringLibrary.R")

# TODO: Do not include segments with lower than 5K bp (see paper)

#
# Load input
#
cd_cnprep()
normal_samples <- load_samples(classes = c("N"), sampleList = "./resources/sampleList.csv")
cytobands <- retrieveCytobands(dir = "./resources/cytoBand.txt")
chromosomeSizes <- readRDS("./resources/chromosomeSizes.rds")
normalSegments <- selectSegmentsWithEvents(events = c("A", "D", "N"), samples = normal_samples, chromosomeSizes = chromosomeSizes, 
                                           dir = "./resources/Project_TUV_12995_B01_SOM_Targeted.2018-03-02/", extension = "cnv.facets.v0.5.2.txt", inSampleFolder = TRUE, 
                                           rescaleInput = TRUE, ampCall = 0.2, delCall = -0.235)
tumor_samples <- load_samples(classes = c("T"), sampleList = "./resources/sampleList.csv")

# Generate norminput argument
norminput <- retrieveNormInput(normalSegments)

for(tumor_samples.i in 1:length(tumor_samples)) {
  sample <- tumor_samples[tumor_samples.i]
  
  print(paste("Analyzing sample", sample))
  
  facets_segment_data <- retrieveFacetsSegments(sample, dir = "./resources/Project_TUV_12995_B01_SOM_Targeted.2018-03-02/")
  facets_snp_data <- retrieveFacetsSnps(sample, dir = "./resources/Project_TUV_12995_B01_SOM_Targeted.2018-03-02/")
  
  # Generate seginput argument
  seginput <- retrieveSegInput(facets_segment_data, sample, cytobands)
  print(paste("Retrieved segment input for sample", sample))
  
  # Generate ratinput argument
  ratinput <- retrieveRatInput(facets_snp_data, sample)
  print(paste("Retrieved ratio input for sample", sample))
  
  # Run CNprep:CNpreprocessing
  segtable <- runCNpreprocessing(seginput = seginput, ratinput = ratinput, norminput = norminput, modelNames = "V") #TODO: Is there a distrib="Grid"?
  print(paste("Produced segtable for sample", sample))
  
  write.table(segtable, paste("./output/", sample, "_segtable.tsv", sep = ""), row.names = F, sep = "\t", quote = FALSE)
  print(paste("Wrote output for sample", sample))
}
