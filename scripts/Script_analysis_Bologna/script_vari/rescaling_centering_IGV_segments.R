# median centering of ascat refitted segments - creation of IGV seg file

library(data.table)
library(tidyverse)

import <- fread("D:/analisi_in_corso/risk_1q_13/IGV/IGV_Rawsegm.seg")

import %>% ggplot(aes(seg.mean)) + geom_density() + xlim(0,5)

samples <- unique( import$ID)
medians <- import %>% group_by(ID) %>% summarise(mediansample= median(seg.mean))


i=7
res <- data.frame()
for ( i in 1:length(samples)) {
  
  sample <- samples[i]
  
  samp_df <- filter(import, ID == sample) 
  
  # samp_df$seg.mean %>% plot
  
  samp_median <- median(samp_df$seg.mean)
  
  samp_df$seg.mean_medAdj <- samp_df$seg.mean - samp_median 
  
  # samp_df$seg.mean_medAdj %>% plot
  
  res <- rbind(res, samp_df)
  
}


# rescale for the scaling factor (ASCAT overestimate)
res$seg.mean_medAdj <- res$seg.mean_medAdj / 1.12


ggplot(res, aes(seg.mean_medAdj)) + geom_density() + xlim(-2,4)

write_tsv(res[,-6], "D:/analisi_in_corso/risk_1q_13/IGV/IGV_Rawsegm_medianCentered_Scaled.seg")

 
# sample <- "MM_cytoscanHD_11162017_880_.CytoScanHD_Array."
# 
# samp_df <- filter(import, ID == sample) 
# samp_df$seg.mean %>% plot
# 
# samp_median <- median(samp_df$seg.mean)
# 
# samp_df$seg.mean_medAdj <- samp_df$seg.mean - samp_median 
# samp_df$seg.mean_medAdj %>% plot
# 



