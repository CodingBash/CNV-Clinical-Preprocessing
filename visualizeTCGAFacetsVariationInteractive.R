#
# This script calculates the SD and RMSE of the TCGA FACETS output from the 
# nonmatching tumor-normal pairs. The statistics were calculate per SNP and
# by segments (sampled at intervals)
#

#
# Import source libraries
#
setwd("~/Git-Projects/Git-Research-Projects/drug-response-prediction")
source("helperFunctions.R")
source("cnProfileVisualizationLibrary.R")
source("visualizeTCGAFacetsVariationLibrary.R")

#
# Visualize SD of SNP for all tumors
#
visualizeSnpSDAllTumor <- function(){
  cd_facets()
  for(tumorId in seq(1, 7)){
    sd_bed <- calculateSnpSD(tumorId)
    visualizeCNProfile(title = paste0("SNP SD Profile for tumorId=", tumorId), facets_snp_data = sd_bed, categories = c("chr14"), save = FALSE, ymin = 0, ymax = 1.5)
  }
}

#
# Visualize RMSE of SNP for all tumors
#
visualizeSnpRMSEAllTumor <- function() {
  cd_facets()
  for(tumorId in seq(1, 7)){
    rmse_bed <- calculateSnpRMSE(tumorId)
    visualizeCNProfile(title = paste0("SNP RMSE Profile for tumorId=", tumorId), facets_snp_data = rmse_bed, categories = c("chr14"), save = FALSE, ymin = 0, ymax = 1.5)
  }
}

#
# Visualize SD of segments for all tumors
#
visualizeSegmentSDAllTumor <- function(){
  cd_facets()
  for(tumorId in seq(1, 7)){
    sd_bed <- calculateSegmentSD(tumorId)
    visualizeCNProfile(title = paste0("Segment SD Profile for tumorId=", tumorId), facets_snp_data = sd_bed, categories = c("chr14"), save = FALSE, ymin = 0, ymax = 1)
  }
}

#
# Visualize RMSE of segments for all tumors
#
visualizeSegmentRMSEAllTumor <- function() {
  cd_facets()
  for(tumorId in seq(1, 7)){
    rmse_bed <- calculateSegmentRMSE(tumorId)
    visualizeCNProfile(title = paste0("Segment RMSE Profile for tumorId=", tumorId), line_data = rmse_bed, categories = c("chr14"), save = FALSE, ymin = 0, ymax = 1)
  }
}

#
# Visualize SD of SNP and segments for all tumors
#
visualizeSegmentAndSnpSDAllTumor <- function() {
  cd_facets()
  for(tumorId in seq(1, 7)){
    seg_sd_bed <- calculateSegmentSD(tumorId)
    snp_sd_bed <- calculateSnpSD(tumorId)
    visualizeCNProfile(title = paste0("Segment SNP SD Profile for tumorId=", tumorId), facets_snp_data = snp_sd_bed, line_data = seg_sd_bed, categories = c("chr14"), save = FALSE, ymin = 0, ymax = 1)
  }
}

#
# Visualize RMSE of SNP and segments for all tumors
#
visualizeSegmentAndSnpRMSEAllTumor <- function() {
  cd_facets()
  for(tumorId in seq(1, 7)){
    seg_rmse_bed <- calculateSegmentRMSE(tumorId)
    snp_rmse_bed <- calculateSnpRMSE(tumorId)
    visualizeCNProfile(title = paste0("Segment SNP RMSE Profile for tumorId=", tumorId), facets_snp_data = snp_rmse_bed, line_data = seg_rmse_bed, categories = c("chr14"), save = FALSE, ymin = 0, ymax = 1)
  }
}

