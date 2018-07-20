source("https://bioconductor.org/biocLite.R")
biocLite("GenomicRanges")

setwd("~/Git-Projects/Git-Research-Projects/drug-response-prediction")
source("helperFunctions.R")
library(GenomicRanges)

samples <- load_samples(classes = c("T"))
sample <- sample[[1]]

cd_doc()
facets_snp_data <- retrieveFacetsSnps(sample, dir = "CSHL/Project_TUV_12995_B01_SOM_Targeted.2018-03-02/")
snp_bed <- snpsToBedFormat(facets_snp_data)
# TODO: First, read in FACETS SNP data into a BED file

# TODO: Second, transform BED file into GRanges object

seqnames <- table(snp_bed[[1]])

gr <- GRanges(
 seqnames = Rle(names(seqnames), as.vector(seqnames)),
 ranges = IRanges(snp_bed[[2]], snp_bed[[3]], names = row.names(snp_bed)),
 score = snp_bed[[4]]
)


# start and end will be the same
print(gr)
