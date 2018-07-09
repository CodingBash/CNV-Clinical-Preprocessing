library(ggplot2)

cd_local <- function() {
  setwd("C:/Users/bbece/Documents/Git-Projects/Git-Research-Projects/drug-response-prediction")
}

cd_local()
samples <- read.table("sampleList.csv", header=T, sep = "\t", stringsAsFactors = F)

# TODO: This is a duplicate from cnvCores.R!!! Modularize
classes <- c("T")
loaded_samples <- c(NA)
loaded_samples.index <- 1
for(sample in samples$Organoids){
  for(class in classes){
    if(substring(sample, 2,2) == class){
      loaded_samples[loaded_samples.index] <- sample
      loaded_samples.index <- loaded_samples.index + 1
      next
    }
  }
}


plots <- list()
for(i in seq(1, length(loaded_samples))){
  sample <- loaded_samples[[i]]
  cd_local()
  
  segtable <- read.table(paste("segClusteringResults/", sample, "_segtable.tsv", sep = ""), sep = "\t", header = TRUE)
  gaussian_comps <- unique(segtable[c("maxzmean", "maxzsigma")])
  
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
  #hist(segtable$segmedian,
  #     main = "Gaussian mixture model of segtable$segmedian,",
  #     col = col,
  #     breaks = tb,
  #     axes = FALSE,
  #     prob = TRUE)
  plt <- ggplot2::ggplot(data = segtable, aes(segtable$segmedian)) + 
    geom_histogram(aes(y = ..density..),
                   breaks = tb,
                   fill = head(col, length(col)-1),
                   alpha = 0.6) + 
    xlim(c(-3,3)) +
    labs(title = "Gaussian mixture model of segtable$segmedian") +
    geom_density(col="black", stat = "density")
  
  
  
  components <- c()
  components.index <- 1
  for(gaussian_comps.index in seq(1, nrow(gaussian_comps))){
    mean <- as.numeric(gaussian_comps[gaussian_comps.index,]$maxzmean)
    sd <- as.numeric(gaussian_comps[gaussian_comps.index,]$maxzsigma)
    a <- rnorm(1000, mean, sd)
    components[[components.index]] <- a
    components.index <- components.index + 1
    plt <- plt + stat_function(fun = dnorm, args = list(mean = mean, sd = sd))
  }
  
  pdf(paste("GMMs/plots_", sample, ".pdf", sep = ""))
  print(plt)
  dev.off()
  #plt <- plt + stat_density(aes(x = components), position = "identity", geom="line")
}

