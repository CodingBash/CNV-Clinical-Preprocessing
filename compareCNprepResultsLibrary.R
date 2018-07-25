#
# Create dataframe of available CNprep runs
#
all_model_specs <- data.frame(stringsAsFactors = FALSE)
cnprep_run <- list()
cnprep_run[[1]] <- data.frame(dir="prev_run1", model="E", minjoin=0.25, ntrial = 10, stringsAsFactors = FALSE)
cnprep_run[[2]] <- data.frame(dir="prev_run3", model="V", minjoin=0.25, ntrial = 10, stringsAsFactors = FALSE)
cnprep_run[[3]] <- data.frame(dir="prev_run_7_15_2018_14", model="E", minjoin=0.25, ntrial = 40, stringsAsFactors = FALSE)
cnprep_run[[4]] <- data.frame(dir="prev_run_7_15_2018_15", model="E", minjoin=0.25, ntrial = 10, stringsAsFactors = FALSE)
cnprep_run[[5]] <- data.frame(dir="prev_run_7_15_2018_16", model="V", minjoin=0.25, ntrial = 10, stringsAsFactors = FALSE)
cnprep_run[[6]] <- data.frame(dir="prev_run_7_15_2018_17", model="E", minjoin=0.01, ntrial = 10, stringsAsFactors = FALSE)
cnprep_run[[7]] <- data.frame(dir="prev_run_7_15_2018_18", model="V", minjoin=0.01, ntrial = 10, stringsAsFactors = FALSE)
cnprep_run[[8]] <- data.frame(dir="prev_run_7_15_2018_19", model="E", minjoin=0.01, ntrial = 50, stringsAsFactors = FALSE)
cnprep_run[[9]] <- data.frame(dir="prev_run_7_15_2018_18", model="V", minjoin=0.01, ntrial = 50, stringsAsFactors = FALSE)
cnprep_run[[10]] <- data.frame(dir="prev_run_7_15_2018_18", model="V", minjoin=0.01, ntrial = 50, stringsAsFactors = FALSE)
cnprep_run[[11]] <- data.frame(dir="prev_run_7_19_2018_1", model="E", minjoin=0.50, ntrial = 40, stringsAsFactors = FALSE)
cnprep_run[[12]] <- data.frame(dir="prev_run_7_19_2018_2", model="E", minjoin=1.00, ntrial = 40, stringsAsFactors = FALSE)
cnprep_run[[13]] <- data.frame(dir="prev_run_7_19_2018_3", model="V", minjoin=1.00, ntrial = 10, stringsAsFactors = FALSE)
all_model_specs <- do.call(rbind, cnprep_run)

#
# Retrieve the segtable from the CNprep::CNpreprocessing output
#
retrieveSegtable <- function(sample, dir = "segClusteringResults/"){
  segtable <- read.table(paste(dir, sample, "_segtable.tsv", sep = ""), sep = "\t", header = TRUE)
  return(segtable)
}

#
# Method to display and compare CNprep results
#
# @param("organoidId") - string of the organoid sample
# @param("bin_start") - start of the x-axis (in bin units)
# @param("bin_end") - end of the x-axis (in bin units)
# @param("model_specs") - list of CNprep runs to compare
# @param("cluster_value") - value of the segtable column with the cluster values
# @param("supplementary_values") - values of the segtable columns with any supplementary data to overlay
# @param("cluster_cols") - color palette for each cluster
# @param("supplementary_cols") - color palette for each supplementary data
# @param("hl") - boolean declaring if plot horizontal lines should appear
#
displayCNprepResults <- function(organoidId, bin_start, bin_end, model_specs, cluster_value = "maxzmean", supplementary_values, cluster_cols, supplementary_cols, hl = TRUE){
  #
  # Set plot parameters
  #
  layout(matrix(seq(1, nrow(model_specs) * 2), nrow(model_specs), 2, byrow = TRUE), 
         widths=c(4,1))
  par(mar=c(2,2,1.25,0))
  
  #
  # Retrieve list of all CNprep segtables (accounting for bin range)
  #
  all_segtables <- lapply(seq(1, nrow(model_specs)), function(model_specs.index, bin_start, bin_end){
    
    segtable <- retrieveSegtable(organoidId, dir = paste0("segClusteringResults/", model_specs[model_specs.index, ]$dir, "/"))
    if(!missing(bin_start) & !missing(bin_end)){
      segtable <- segtable[segtable$start >= bin_start & segtable$end <= bin_end,]  
    } else if (!missing(bin_start)){
      segtable <- segtable[segtable$start >= bin_start,]  
    } else if (!missing(bin_end)) {
      segtable <- segtable[segtable$end <= bin_end,]  
    }
    return(segtable)
  }, bin_start, bin_end)
  
  #
  # Generate dataframe with all CNprep segtables
  #
  binded_segtables <- do.call(rbind, all_segtables)

  #
  # Determine plot ranges
  #
  values <- list()
  if(!missing(cluster_value)){
    values <- c(values, cluster_value)
  }
  if(!missing(supplementary_values)){
    values <- c(values, supplementary_values)
  }
  ymargin <- 0.1
  yplot <- 0.5
  yvalues <- sapply(values, function(value){binded_segtables[[value]]})
  xrange <- range(binded_segtables$start,binded_segtables$end)
  yrange <- range(yvalues - ymargin, yvalues + ymargin + 0.5)
    
  #
  # Iterate through each CNprep segtable
  #
  for(segtable.i in seq_along(all_segtables)){
    #
    # Store necessary data information
    #
    segtable <- all_segtables[[segtable.i]]
    model_spec <- model_specs[segtable.i, ]
    
    #
    # Set plot information
    #
    as.list(model_spec)
    title <- paste0("dir=", model_spec$dir, " model=", model_spec$model, " minjoin=", model_spec$minjoin, " ntrial=", model_spec$ntrial)
    plot(xrange, yrange, main = title, type="n", xlab = "", ylab = "")
    legend_values <- c()
    legend_col <- c()
    
    #
    # Display supplementary value segments
    #
    if(!missing(supplementary_values)){
      for(supplementary_value.index in seq_along(supplementary_values)){
        segments(x0 = segtable$start, x1 = segtable$end, y0 = segtable[[supplementary_values[[supplementary_value.index]]]], y1 = segtable[[supplementary_values[[supplementary_value.index]]]],
                 col = supplementary_cols[[supplementary_value.index]], lty = par("lty"), lwd = par("lwd"))
      }
      
      legend_values <- c(legend_values, supplementary_values)
      legend_col <- c(legend_col, head(supplementary_cols, length(supplementary_values)))
    }
    
    #
    # Display clusters
    #
    if(!missing(cluster_value)){
        clusters <- split(segtable, f=segtable[[cluster_value]])
        for(cluster.index in seq_along(clusters)){
          segments(x0 = clusters[[cluster.index]]$start, x1 = clusters[[cluster.index]]$end, y0 = clusters[[cluster.index]][[cluster_value]], y1 = clusters[[cluster.index]][[cluster_value]],
                   col = cluster_cols[[cluster.index]], lty = par("lty"), lwd = par("lwd"))
        }
        legend_values <- c(legend_values, paste0(cluster_value,"#", unlist(lapply(names(clusters), function(cluster){ return(substr(cluster, 1, 5))}))))
        legend_col <- c(legend_col, head(cluster_cols, length(names(clusters))))
    }
    
    #
    # Display horizontal lines
    #
    if(hl == TRUE){
      abline(h=0)
      abline(h=0.5, col = "#4F4F4F")
      abline(h=-0.5, col = "#4F4F4F")
      abline(h=1, col = "#BABABA")
      abline(h=-1, col = "#BABABA")
    }
    
    #
    # Set legend
    #
    plot.new()
    legend("left", legend=legend_values,
           col = legend_col, seg.len = 0.5, lty=1, cex=1, pt.cex = 1, xpd = TRUE, y.intersp=1, x.intersp=.1)
  }
}