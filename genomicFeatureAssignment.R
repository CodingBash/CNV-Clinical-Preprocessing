#source("https://bioconductor.org/biocLite.R")
#biocLite("GenomicRanges")

source("helperFunctions.R")
library(GenomicRanges)
library(ggplot2) 
library(reshape) 
library(reshape2)
library(made4)
library(cluster)
library(spatstat) # "im" function 

#
# Load sample to retrieve feature set for
#

samples <- load_samples(classes = c("T", "M"))

retrieveCores <- function(dir){
  return(read.table(dir, header = FALSE, sep = "\t", stringsAsFactors = FALSE))
}

retrieveTrainingSet <- function(loaded_samples, Acores, Dcores, binDir = "CSHL/Project_TUV_12995_B01_SOM_Targeted.2018-03-02/", removeNAs = TRUE){
  training_set <- data.frame(stringsAsFactors = FALSE)  
  for(sample in samples){
    #
    # Retrieve necessary data for feature value calculation (bins with cnlr and COREs)
    #
    facets_snp_data <- retrieveFacetsSnps(sample, dir = binDir)
    snp_bed <- snpsToBedFormat(facets_snp_data)
    
    #
    # Preprocess COREs
    #
    seqnames_A <- table(Acores[[1]])
    gr_A <- GRanges(
      seqnames = Rle(names(seqnames_A), as.vector(seqnames_A)),
      ranges = IRanges(Acores[[2]], Acores[[3]], names = row.names(Acores)),
      event=rep("A", nrow(Acores))
    )
    seqnames_D <- table(Dcores[[1]])
    gr_D <- GRanges(
      seqnames = Rle(names(seqnames_D), as.vector(seqnames_D)),
      ranges = IRanges(Dcores[[2]], Dcores[[3]], names = row.names(Dcores)),
      event=rep("D", nrow(Dcores))
    )
    # TODO: Do operations on GRanges object to simplify, if need be (i.e. reduction)
    gr <- c(gr_A, gr_D)
    coreDf <- as.data.frame(gr, row.names = seq_along(gr$event))
    coreDf$cnlr <- NA
    
    #
    # Determine value of each CORE feature (cnlr)
    #
    coreDf[, 7] <- sapply(seq(1, nrow(coreDf)), function(core.index){
      # Get snps with same chromosome as core, and in between the core region
      snps <- snp_bed[snp_bed[[1]] == coreDf[core.index, 1] & snp_bed[[2]] >= coreDf[core.index, 2] & snp_bed[[3]] <= coreDf[core.index, 3], ]
      # Calculate and assign median
      return(median(snps$value))
    } )
    
    # Convert back to Granges object - this is our final feature object
    featureSet <- makeGRangesFromDataFrame(coreDf, keep.extra.columns = TRUE)
    
    sampleTrainingSet <- data.frame(cnlr = coreDf[,7])
    sampleTrainingSet$coreId <- rownames(coreDf)
    sampleTrainingSet$sampleId <- sample
    
    training_set <- rbind(training_set, sampleTrainingSet)
  }
  if(removeNAs == TRUE){
    training_set <- training_set[-which(is.na(training_set$cnlr)),]
  }
  return(training_set)
}

visualizeUnclusteredHeatmap <- function(training_set){
  ggplot(data = training_set, aes(x = coreId, y = sampleId)) + 
    geom_tile(aes(fill = cnlr), color = "white", size = 1) + 
    scale_fill_gradient2(low = "blue", mid="white", high = "tomato") + 
    xlab("core ID") + 
    theme_grey(base_size = 10) + 
    ggtitle("Heatmap (ggplot)") + 
    theme(axis.ticks = element_blank(), 
          panel.background = element_blank(), 
          plot.title = element_text(size = 12, colour = "gray50")) 
}

clusterTrainingSet <- function(training_set, visualize = FALSE){
  # Unmelt training set for correlation analysis
  training_set_matrix <- dcast(data = training_set,formula = sampleId~coreId,fun.aggregate = sum,value.var = "cnlr")
  sampleIds <- training_set_matrix$sampleId
  training_set_matrix <- training_set_matrix[,-c(1)]
  training_set_matrix <- t(training_set_matrix)
  colnames(training_set_matrix) <- sampleIds
 
  # Calculate distance matrix
  corRaw <- cor(training_set_matrix)
  dissimilarity <- 1 - corRaw
  distance.sample <- as.dist(dissimilarity)
  distance.core <- as.dist(t(dissimilarity))
  
  # Run hierarchical clustering
  hc.sample <- hclust(distance.sample)
  hc.core <- hclust(distance.core)
  
  if(visualize == TRUE){
    color.palette  <- colorRampPalette(c("blue", "white", "tomato"))(n=600)
    heatmap.2(t(training_set_matrix),  trace="none", dendrogram="row", density.info = 'none', scale='none', col = color.palette)
  }
  
  return(hc.sample)  
}