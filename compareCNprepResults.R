setwd("~/Git-Projects/Git-Research-Projects/drug-response-prediction")
source("helperFunctions.R")
library(RColorBrewer)

all_model_specs <- data.frame(stringsAsFactors = FALSE)
cnprep_run_1 <- data.frame(dir="prev_run1", model="E", minjoin=0.25, ntrial = 10, stringsAsFactors = FALSE)
cnprep_run_2 <- data.frame(dir="prev_run3", model="V", minjoin=0.25, ntrial = 10, stringsAsFactors = FALSE)
cnprep_run_3 <- data.frame(dir="prev_run_7_19_2018_1", model="E", minjoin=0.50, ntrial = 40, stringsAsFactors = FALSE)
cnprep_run_4 <- data.frame(dir="prev_run_7_19_2018_2", model="E", minjoin=1.00, ntrial = 40, stringsAsFactors = FALSE)
cnprep_run_5 <- data.frame(dir="prev_run_7_19_2018_3", model="V", minjoin=1.00, ntrial = 10, stringsAsFactors = FALSE)
all_model_specs <- rbind(cnprep_run_1, cnprep_run_2, cnprep_run_3, cnprep_run_4, cnprep_run_5)

displayCNResults <- function(organoidId, bin_start, bin_end, model_specs, values, cols, hl){
  par(mfrow=c(nrow(model_specs),1)) 
  par(mar=c(2,2,1.25,1))
  #par(cex.main = .5, cex.axis = 0.5)
  
  #
  # Calculate plot range
  #
  all_segtables <- lapply(seq(1, nrow(model_specs)), function(model_specs.index){
    segtable <- retrieveSegtable(organoidId, dir = paste0("segClusteringResults/", model_specs[model_specs.index, ]$dir, "/"))
    segtable <- segtable[segtable$start >= bin_start & segtable$end <= bin_end,]
    return(segtable)
  })
  
  binded_segtables <- do.call(rbind, all_segtables)
  
  ymargin <- 0.1
  yvalues <- sapply(values, function(value){binded_segtables[[value]]})
  xrange <- range(binded_segtables$start,binded_segtables$end)
  yrange <- range(yvalues - ymargin, yvalues + ymargin)
    
  for(segtable.i in seq_along(all_segtables)){
    segtable <- all_segtables[[segtable.i]]
    model_spec <- model_specs[segtable.i, ]
    as.list(model_spec)
    title <- paste0("dir=", model_spec$dir, " model=", model_spec$model, " minjoin=", model_spec$minjoin, " ntrial=", model_spec$ntrial)
    plot(xrange, yrange, main = title, type="n", xlab = "", ylab = "")
    for(value.i in seq_along(values)){
      segments(x0 = segtable$start, x1 = segtable$end, y0 = segtable[[values[[value.i]]]], y1 = segtable[[values[[value.i]]]],
               col = cols[[value.i]], lty = par("lty"), lwd = par("lwd"))
    }
    if(hl == TRUE){
      abline(h=0)
      abline(h=0.5, col = "#4F4F4F")
      abline(h=-0.5, col = "#4F4F4F")
      abline(h=1, col = "#BABABA")
      abline(h=-1, col = "#BABABA")
    }
  }
}

cluster_cols <- brewer.pal(n = 7, name="Set1")
supplementary_cols <- brewer.pal(n = 7, name="Set3")
displayCNResults(organoidId= "hT1", bin_start = 20000, bin_end = 40000, model_specs = rbind(cnprep_run_4,cnprep_run_2 ), cluster_value = "maxzmean", supplementary_values = c("segmedian", "mediandev"), cluster_cols = cluster_cols, supplementary_cols =  supplementary_cols, hl = TRUE)
