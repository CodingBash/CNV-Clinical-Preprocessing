#
# Save seginputs from samples in the format accepted by CNprep::CNpreprocessing()
#

#
# Import functions
#
setwd("~/Git-Projects/Git-Research-Projects/drug-response-prediction")
source("helperFunctions.R")
source("segmentClusteringLibrary.R")

#
# Load samples to save
#
cd_local()
input_samples <- load_samples(classes = c("T"), sampleList = "sampleList.csv")

#
# Iterate through samples
#
for(input_sample in input_samples){
  #
  # Generate seginput
  #
  cd_doc()
  facets_segment_data <- retrieveFacetsSegments(input_sample, dir = "CSHL/Project_TUV_12995_B01_SOM_Targeted.2018-03-02/")
  seginput <- retrieveSegInput(facets_segment_data, cytobands)
  
  #
  # Write to directory
  #
  cd_local()
  write.table(seginput, paste("./segInputs/", input_sample, "_seginput.tsv", sep = ""), row.names = F, sep = "\t", quote = FALSE)
  print(paste("Wrote output for sample", input_sample))
  
}
