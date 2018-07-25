
#
# Load sources
#
setwd("~/Documents/Git-Projects/Git-Research-Projects/drug-response-prediction/") # Set working directory to where the scripts are located at
source("compareCNprepResultsLibrary.R") # Import visualization library
library(RColorBrewer) # Import brewer for coloring
library(repr)

cluster_cols <- brewer.pal(n = 7, name="Dark2") # Set colors for clusters. Let n > 5
supplementary_cols <- brewer.pal(n = 7, name="Set2") # Set colors for suppl. Let n >= number of suppl values

options(repr.plot.width=15, repr.plot.height=7)


View(all_model_specs)

print(all_model_specs)

displayCNprepResults(organoidId= "hT1", model_specs = all_model_specs[c(1,12), ], cluster_value = "maxzmean", cluster_cols = cluster_cols) # Display

displayCNprepResults(organoidId= "hT1", model_specs = all_model_specs[c(1,12), ], cluster_value = "maxzmean", supplementary_values = c("mediandev"), cluster_cols = cluster_cols, supplementary_cols =  supplementary_cols)


displayCNprepResults(organoidId= "hT1", bin_start = 20000, bin_end = 40000, model_specs = all_model_specs[c(1,12), ], cluster_value = "maxzmean", supplementary_values = c("mediandev"), cluster_cols = cluster_cols, supplementary_cols =  supplementary_cols)

displayCNprepResults(organoidId= "hT1", bin_start = 20000, bin_end = 40000, model_specs = all_model_specs[c(1,12), ], cluster_value = "maxzmean", supplementary_values = c("mediandev"), cluster_cols = cluster_cols, supplementary_cols =  supplementary_cols, hl = FALSE)
