#
# Load sources
#
setwd("~/Git-Projects/Git-Research-Projects/drug-response-prediction/") # Set working directory to where the scripts are located at
source("compareCNprepResultsLibrary.R") # Import visualization library
library(RColorBrewer) # Import brewer for coloring

#
# Set colors
#
cluster_cols <- brewer.pal(n = 7, name="Dark2") # Set colors for clusters. Let n > 5
supplementary_cols <- brewer.pal(n = 7, name="Set2") # Set colors for suppl. Let n >= number of suppl values


# Example 1: Display CNprep results for organoid "hT1" for CNprep runs #1,12. 
displayCNprepResults(organoidId= "hT1", model_specs = all_model_specs[c(1,12), ], cluster_value = "maxzmean", cluster_cols = cluster_cols) # Display

# Example 2: Extended off of example 1 by also viewing the mediandev of segments too
displayCNprepResults(organoidId= "hT1", model_specs = all_model_specs[c(1,12), ], cluster_value = "maxzmean", supplementary_values = c("mediandev"), cluster_cols = cluster_cols, supplementary_cols =  supplementary_cols)

# Example 3: Extended off of example 2 by also viewing all segments between bins 20000 and 40000
displayCNprepResults(organoidId= "hT1", bin_start = 20000, bin_end = 40000, model_specs = all_model_specs[c(1,12), ], cluster_value = "maxzmean", supplementary_values = c("mediandev"), cluster_cols = cluster_cols, supplementary_cols =  supplementary_cols)