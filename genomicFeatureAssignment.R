#source("https://bioconductor.org/biocLite.R")
#biocLite("GenomicRanges")

setwd("~/Git-Projects/Git-Research-Projects/drug-response-prediction")
source("helperFunctions.R")
library(GenomicRanges)
library(ggplot2) 
library(reshape) 
library(reshape2)
library(cluster)

#
# Load sample to retrieve feature set for
#
training_set <- data.frame(stringsAsFactors = FALSE)
samples <- load_samples(classes = c("T", "M"))
for(sample in samples){
  
  #
  # Retrieve feature labels (COREs) analgous to all samples
  cd_local("hT_results/finalCores/prev_run_7_3_2018")
  Acores <- read.table(paste0("AfinalCoresBP_subtractFraction.bed"), header = FALSE, sep = "\t", stringsAsFactors = FALSE) 
  Dcores <- read.table(paste0("DfinalCoresBP_subtractFraction.bed"), header = FALSE, sep = "\t", stringsAsFactors = FALSE) 
  
  #
  # Retrieve necessary data for feature value calculation (bins with cnlr and COREs)
  #
  cd_doc()
  facets_snp_data <- retrieveFacetsSnps(sample, dir = "CSHL/Project_TUV_12995_B01_SOM_Targeted.2018-03-02/")
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
  print(featureSet)
  
  sampleTrainingSet <- data.frame(cnlr = coreDf[,7])
  sampleTrainingSet$coreId <- rownames(coreDf)
  sampleTrainingSet$sampleId <- sample
  
  training_set <- rbind(training_set, sampleTrainingSet)
 
}

training_set <- training_set[-which(is.na(training_set$cnlr)),]
ggplot(data = training_set, aes(x = coreId, y = sampleId)) + 
  geom_tile(aes(fill = cnlr), color = "white", size = 1) + 
  scale_fill_gradient2(low = "blue", mid="white", high = "tomato") + 
  xlab("core ID") + 
  theme_grey(base_size = 10) + 
  ggtitle("Heatmap (ggplot)") + 
  theme(axis.ticks = element_blank(), 
        panel.background = element_blank(), 
        plot.title = element_text(size = 12, colour = "gray50")) 

# Unmelt training set for correlation analysis
training_set_matrix <- dcast(data = training_set,formula = sampleId~coreId,fun.aggregate = sum,value.var = "cnlr")
rownames <- training_set_matrix$sampleId
training_set_matrix <- training_set_matrix[,-c(1)]
training_set_matrix <- t(training_set_matrix)
colnames(training_set_matrix) <- rownames
corRaw <- cor(training_set_matrix)

library(spatstat) # "im" function 
plot(im(corRaw[nrow(corRaw):1,]), main="Correlation Matrix Map")

dissimilarity <- 1 - corRaw
distance.sample <- as.dist(dissimilarity)
distance.core <- as.dist(t(dissimilarity))

hc.sample <- hclust(distance.sample)
hc.core <- hclust(distance.core)

# draw heatmap for first cluster
plot(hc.sample)
color.palette  <- colorRampPalette(c("blue", "white", "tomato"))(n=600)
heatmap.2(t(training_set_matrix),  trace="none", dendrogram="row", density.info = 'none', scale='none', col = color.palette)

plot(hclust(distance), 
     main="Dissimilarity = 1 - Correlation", xlab="")
