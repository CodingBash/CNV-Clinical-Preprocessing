install.packages("ggplot2")
library(ggplot2)
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
    if(tb[tb.index] >= min(segmedians) - (interval) & tb[tb.index] <= max(segmedians) + (interval)){
      col[tb.index] <- color  
    } 
  }
  color <- color + 1
}

col[which(col==1)] <- "black"
col[which(col==2)] <- "firebrick1"

col[which(col==3)] <- "gold"
col[which(col==4)] <- "purple"
col[which(col==5)] <- "darkolivegreen1"
hist(segtable$segmedian,
     main = "Gaussian mixture model of segtable$segmedian,",
     col = col,
     breaks = tb,
     axes = FALSE,
     prob = TRUE)
ggplot(data = segtable, aes(segtable$segmedian)) + 
  geom_histogram(aes(y = ..density..),
                 breaks = tb,
                 fill = head(col, length(tb)-1),
                 alpha = 0.6) + 
  labs(title = "Gaussian mixture model of segtable$segmedian") +
  geom_density(col="black")


axis(side = 1, at=tb)
axis(side = 2, at=seq(0, 20, 1))
lines(density(segtable$segmedian))
