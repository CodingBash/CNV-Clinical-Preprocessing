#
# This script is a set of functions used to visualize a GMM from the CNprep::CNpreprocessing result
#


#
# Load visualization libraries
#
library(ggplot2)
library(reshape2)

#
# Retrieve the segtable from the CNprep::CNpreprocessing output
#
retrieveSegtable <- function(sample, dir = "segClusteringResults/"){
  segtable <- read.table(paste(dir, sample, "_segtable.tsv", sep = ""), sep = "\t", header = TRUE)
  return(segtable)
}

#
# Retrieve the seginput from the CNprep::CNpreprocessing input
#
retrieveSeginput <- function(sample, dir = "segInput/"){
  seginput <- read.table(paste(dir, sample, "_seginput.tsv", sep = ""), sep = "\t", header = TRUE)
  return(seginput)
}

#
# Display the GMM visualization using the CNprep::CNpreprocessing output segtable
# TODO: Save location is hardcoded
#
displayGMM <- function(segtable, sample, print = FALSE, save = FALSE){
  
  # Retrieve dataframe of GMM components (w/ mean and sigma)
  gaussian_comps <- unique(segtable[c("maxzmean", "maxzsigma")])
  
  #
  # Set histogram values
  #
  interval <- 0.1 # Bin interval length
  min <- -2.5 # Bin minimum value
  max <- -min # Bin maximum value
  tb <- seq(min, max, interval) # List of all bins
  colors <- c("black", "red", "orange", "yellow", "blue", "purple", "gray") # Colors of histogram components
  
  #
  # Retrieve segment$mediandev for each cluster in separate lists
  #
  cluster_list <- list()
  for(gaussian_comps.index in seq(1, nrow(gaussian_comps))){
    cluster_medians <- list(segtable[gaussian_comps[gaussian_comps.index, ]$maxzmean == segtable$maxzmean & gaussian_comps[gaussian_comps.index, ]$maxzsigma == segtable$maxzsigma, ]$mediandev)
    cluster_list[gaussian_comps.index] <- cluster_medians
  }
  
  # Melt the 2D cluster list into an appropriate format for the ggplot2::geom_histogram function
  hist_data <- melt(cluster_list)
  
  #
  # Create the plot
  # TODO: geom_density is not working - currently creates density within the scope of each GMM component instead of whole model
  #
  plt <- ggplot(data = hist_data, aes(x=value, fill = as.factor(L1))) + 
    geom_histogram(aes(y = ..density..),
                   breaks = tb,
                   position = "stack",
                   alpha = 0.6) +
#   geom_density(data = hist_data, aes(x=value), col="black", stat = "density", position = "identity") + 
    scale_fill_manual(values=tail(colors, length(colors) - 1)) + 
    xlim(c(min,max)) +
    labs(title = paste("Gaussian mixture model of segtable$mediandev for sample", sample))
  
  # Add the GMM component functions to the plot
  for(gaussian_comps.index in seq(1, nrow(gaussian_comps))){
    mean <- as.numeric(gaussian_comps[gaussian_comps.index,]$maxzmean)
    sd <- as.numeric(gaussian_comps[gaussian_comps.index,]$maxzsigma)
    plt <- plt + stat_function(fun = dnorm, n=1000, args = list(mean = mean, sd = sd), col = colors[gaussian_comps.index + 1])
  }
  
  #
  # Print and/or save the plot
  #
  if(print == TRUE){
    print(plt)
  }
  if(save == TRUE){
    ggsave(filename=paste("GMMs/mediandev/plots_", sample, ".pdf", sep = ""), plot=plt, width=16, height=9, units="in")
  }
}

#
# Display the segtable histogram
# TODO: Save location is hardcoded
#
displaySegtableHistogram <- function(segtable, sample, print = FALSE, save = FALSE){
  #
  # Set histogram values
  #
  interval <- 0.1 # Bin interval length
  min <- -2.5 # Bin minimum value
  max <- -min # Bin maximum value
  tb <- seq(min, max, interval) # List of all bins
  print(paste0("Display segtable histogram for tumor=", sample, " which contains the segment_count=", nrow(segtable)))  
  plt <- ggplot(data = segtable, aes(segtable$mediandev)) + 
    geom_histogram(breaks = tb,
                   alpha = 0.6) + 
    labs(title = "Gaussian mixture model of segtable$mediandev")
  
  #
  # Print and/or save the plot
  #
  if(print == TRUE){
    print(plt)
  }
  if(save == TRUE){
    #ggsave(filename=paste("GMMs/mediandev/plots_", sample, ".pdf", sep = ""), plot=plt, width=16, height=9, units="in")
    # TODO: saving
  }
}

#
# Display the seginput histogram
# TODO: Save location is hardcoded
#
displaySeginputHistogram <- function(seginput, sample, print = FALSE, save = FALSE){
  
  #
  # Set histogram values
  #
  interval <- 0.1 # Bin interval length
  min <- -2.5 # Bin minimum value
  max <- -min # Bin maximum value
  tb <- seq(min, max, interval) # List of all bins
  print(paste0("Display seginput histogram for tumor=", sample, " which contains the segment_count=", nrow(seginput)))  
  plt <- ggplot(data = seginput, aes(seginput$seg.median)) + 
    geom_histogram(breaks = tb,
                   alpha = 0.6) + 
    labs(title = "Gaussian mixture model of segtable$mediandev")
  
  #
  # Print and/or save the plot
  #
  if(print == TRUE){
    print(plt)
  }
  if(save == TRUE){
    #ggsave(filename=paste("GMMs/mediandev/plots_", sample, ".pdf", sep = ""), plot=plt, width=16, height=9, units="in")
    # TODO: Saving
  }
}
