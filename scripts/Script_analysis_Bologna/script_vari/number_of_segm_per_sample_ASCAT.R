# inspect the number of segments per sample from refitted ASCAT files!

library(data.table)
library(tidyverse)

setwd("D:/analisi_in_corso/risk_1q_13/ASCAT/Refitted_ASCAT_segm")

segmHD1 <- fread("refitted_ASCAT_segments_batchHD1.tsv")
segmHD2 <- fread("refitted_ASCAT_segments_batchHD2.tsv")
segmHDF <- fread("refitted_ASCAT_segments_batchHDF.tsv")
segm6.0 <- fread("refitted_ASCAT_segments_batchSNP6.tsv")

par(mfrow=c(4,1))

HD1 <- segmHD1$sample %>% table
plot(HD1, main= "HD 1")

HD2 <- segmHD2$sample %>% table
plot(HD2, main= "HD 2")

SNP6 <- segm6.0$sample %>% table
plot(SNP6, main= "SNP6")

HDF <- segmHDF$sample %>% table
plot(HDF, main= "HD F")



