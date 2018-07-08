sample <- "hT1"

segtable <- read.table(paste("segClusteringResults/", sample, "_segtable.tsv", sep = ""), sep = "\t", header = TRUE)

interval <- 0.1
tb <- seq(-3, 3, interval)
col <- rep(1, length(tb))
means <- sort(unique(segtable$maxzmean))
color <- 2
for(mean in means){
  segmedians <- segtable[segtable$maxzmean == mean, ]$segmedian
  print(paste(mean, "cluster between", min(segmedians), "and", max(segmedians)))
  for(tb.index in seq(1, length(tb))){
    if(tb[tb.index] >= min(segmedians) - (interval/2) & tb[tb.index] <= max(segmedians) + (interval/2)){
      col[tb.index] <- color  
    } 
  }
  color <- color + 1
}

col[which(col==1)] <- "black"
col[which(col==2)] <- "firebrick1"

col[which(col==3)] <- "gold"
col[which(col==4)] <- "black"
col[which(col==5)] <- "darkolivegreen1"
hist(segtable$segmedian,
     main = "Gaussian mixture model of segtable$segmedian,",
     col = col,
     breaks = tb,
     axes = FALSE)
axis(side = 1, at=tb)
axis(side = 2, at=seq(0, 20, 1))
#lines(density(segtable$segmedian))
