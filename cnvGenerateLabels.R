#
# This script iterates through the directory will all sample files and creates a CSV containing the names of all samples
# to support downstream analysis scripts.
# Output file is named "sampleList.csv"
#

source("helperFunctions.R")

cd_doc()
files <- list.dirs(path = "CSHL/Project_TUV_12995_B01_SOM_Targeted.2018-03-02/", full.names = FALSE, recursive = FALSE)

sample_list <- c()
sample_list.index <- 1

str(sample_list)
for(file in files){
  sample_text <- substring(file, 1,6)
  if(toupper(sample_text) == "SAMPLE"){
   sample_list[sample_list.index] <- substring(file, 8, nchar(file)) 
   sample_list.index <- sample_list.index + 1
  }
}

print(sample_list)


sample_df <- data.frame(
  Organoids = c(sample_list),
  stringsAsFactors = FALSE
)

print(sample_df$Organoids)

cd_local()
write.csv(sample_df, "sampleList.csv", row.names = FALSE)