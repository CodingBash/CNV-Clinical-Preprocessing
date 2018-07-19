
#
# Import source libraries
#
setwd("~/Documents/Git-Projects/Git-Research-Projects/drug-response-prediction")
source("helperFunctions.R")
source("cnProfileVisualizationLibrary.R")
source("visualizeTCGAFacetsVariationLibrary.R")



visualizeSegmentAndSnpSDAllTumor <- function() {
  setwd(paste0("~/Documents/Git-Projects/Git-Research-Projects/FACETS_nonmatching_test/"))
  for(tumorId in seq(2, 2)){
    seg_sd_bed <- calculateSegmentSD(tumorId)
    snp_sd_bed <- calculateSnpSD(tumorId)
    visualizeCNProfile(title = paste0("Segment SNP SD Profile for tumorId=", tumorId), facets_snp_data = snp_sd_bed, line_data = seg_sd_bed, save = FALSE, ymin = 0, ymax = 1)
  }
}
#facets_snp_data = snp_sd_bed, 
visualizeSegmentAndSnpSDAllTumor()

visualizeSegmentAndSnpRMSEAllTumor <- function() {
  setwd(paste0("~/Documents/Git-Projects/Git-Research-Projects/FACETS_nonmatching_test/"))
  for(tumorId in seq(1, 7)){
    seg_rmse_bed <- calculateSegmentRMSE(tumorId)
    snp_rmse_bed <- calculateSnpRMSE(tumorId)
    visualizeCNProfile(paste0("Segment SNP RMSE Profile for tumorId=", tumorId), facets_snp_data = snp_rmse_bed, line_data = seg_rmse_bed, categories = c("chr14"), save = FALSE, ymin = 0, ymax = 1)
  }
}

visualizeSegmentAndSnpSDAllTumor()



