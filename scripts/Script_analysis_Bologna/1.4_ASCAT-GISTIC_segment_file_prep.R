# from ASCAT to GISTIC format script

# ====== data import ======
library(data.table)
library(tidyverse)

setwd("D:/analisi_in_corso/risk_1q_13/ASCAT/Refitted_ASCAT_segm")

segmHD1 <- fread("refitted_ASCAT_segments_batchHD1.tsv")
segmHD2 <- fread("refitted_ASCAT_segments_batchHD2.tsv")
segm6.0 <- fread("refitted_ASCAT_segments_batchSNP6.tsv")
# segmHDF <- fread("refitted_ASCAT_segments_batchHDF.tsv")



#======= Counting of number of probes per segment (required for GISTIC format) ===========

# import PROBES INFO for both types of arrays (pre-generated from rawcopy)
probesHD <- fread("D:/analisi_in_corso/risk_1q_13/info_files/probePos_HD.txt")
probes6.0 <- fread("D:/analisi_in_corso/risk_1q_13/info_files/probePos_SNP6.txt")


# HD1
num_probes <- c()
for (i in 1:nrow(segmHD1)){
  if (i==1){ 
    t0=Sys.time()
  }
  if(i%%100==0) {
    tx <- Sys.time()
    print(paste0(i, " - elapsed time: ", tx- t0))
  }
  s <- segmHD1[i,]
  p <- probesHD %>% filter(chrs==s$chr, 
                           pos >= s$startpos, 
                           pos < s$endpos )
  n <- nrow(p)
  num_probes <- append(num_probes, n)
}

segmHD1$num_probes <- num_probes
num_probes %>% plot()

# HD2
num_probesHD2 <- c()
for (i in 1:nrow(segmHD2)){
  if (i==1){ 
    t0=Sys.time()
  }
  if(i%%100==0) {
    tx <- Sys.time()
    print(paste0(i, " - elapsed time: ", tx- t0))
  }
  s <- segmHD2[i,]
  p <- probesHD %>% filter(chrs==s$chr, 
                           pos >= s$startpos, 
                           pos < s$endpos )
  n <- nrow(p)
  num_probesHD2 <- append(num_probesHD2, n)
}

segmHD2$num_probes <- num_probesHD2


# # HDF
# num_probesHDF <- c()
# for (i in 1:nrow(segmHDF)){
#   if (i==1){ 
#     t0=Sys.time()
#   }
#   if(i%%100==0) {
#     tx <- Sys.time()
#     print(paste0(i, " - elapsed time: ", tx- t0))
#   }
#   s <- segmHDF[i,]
#   p <- probesHD %>% filter(chrs==s$chr, 
#                            pos >= s$startpos, 
#                            pos < s$endpos )
#   n <- nrow(p)
#   num_probesHDF <- append(num_probesHDF, n)
# }
# 
# segmHDF$num_probes <- num_probesHDF


# SNP 6.0
num_probesSNP6 <- c()
for (i in 1:nrow(segm6.0)){
  if (i==1){ 
    t0=Sys.time()
  }
  if(i%%100==0) {
    tx <- Sys.time()
    print(paste0(i, " - elapsed time: ", tx- t0))
  }
  s <- segm6.0[i,]
  p <- probes6.0 %>% filter(chrs==s$chr, 
                           pos >= s$startpos, 
                           pos < s$endpos )
  n <- nrow(p)
  num_probesSNP6 <- append(num_probesSNP6, n)
}

segm6.0$num_probes <- num_probesSNP6


#================== merging the 4 batches segments ====================

segm <- rbind(segmHD1, 
              segmHD2, 
              # segmHDF, 
              segm6.0)

#======== integer CN to log2R conversion (required for GISTIC format) ==========

segm$CNvalue_refit_logR <- log2(segm$CNvalue_refit) - 1

segm$CNvalueRaw_refit_logR <- log2(segm$CNvalueRaw_refit) - 1

segm %>% ggplot(aes(CNvalue_refit_logR)) + geom_density() + xlim(c(-2,2)) + geom_vline(xintercept = 0.585)

segm %>% ggplot(aes(CNvalueRaw_refit_logR)) + geom_density() + xlim(c(-2,2)) + geom_vline(xintercept = 0.585)


#================ BONUS: IGV files creation ====================

# LogR Raw copy number (refitted for WGD)
IGV_Rawsegm_logR <- data.frame(ID = segm$sample, chrom =segm$chr, loc.start= segm$startpos,
                                loc.end = segm$endpos, num.mark	= segm$num_probes,
                                seg.mean = segm$CNvalueRaw_refit_logR)
readr::write_tsv(IGV_Rawsegm_logR, "C:/Users/mm_gr/Desktop/IGV_Rawsegm_logR.seg" )

# LogR copy number (refitted for WGD)
IGV_Segm_logR <- data.frame(ID = segm$sample, chrom =segm$chr, loc.start= segm$startpos,
                             loc.end = segm$endpos, num.mark	= segm$num_probes,
                             seg.mean = segm$CNvalue_refit_logR)
readr::write_tsv(IGV_Segm_logR, "C:/Users/mm_gr/Desktop/IGV_Segm_logR.seg" )

# Integer Raw copy number (refitted for WGD)
IGV_Rawsegm <- data.frame(ID = segm$sample, chrom =segm$chr, loc.start= segm$startpos,
                               loc.end = segm$endpos, num.mark	= segm$num_probes,
                               seg.mean = segm$CNvalueRaw_refit)
readr::write_tsv(IGV_Rawsegm, "C:/Users/mm_gr/Desktop/IGV_Rawsegm.seg" )

# Integer copy number (refitted for WGD)
IGV_Segm <- data.frame(ID = segm$sample, chrom =segm$chr, loc.start= segm$startpos,
                            loc.end = segm$endpos, num.mark	= segm$num_probes,
                            seg.mean = segm$CNvalue_refit)
readr::write_tsv(IGV_Segm, "C:/Users/mm_gr/Desktop/IGV_Segm.seg" )


#================= GISTIC input file creation ======================

#______ create segm_raw file ____________ MAIN

GISTIC.segm<- data.frame( 
                       "Sample" = segm$sample,
                       "Chromosome" = segm$chr ,
                       "Start Position" = segm$startpos, #+1, # if necessarry manually add 1 base to avoid segment overlapping (GISTIC error)
                       "End Position" =segm$endpos ,
                       "Num markers" = segm$num_probes ,
                       "Seg.CN"= segm$CNvalueRaw_refit_logR
                       )

GISTIC.segm.no_XYM<- filter(GISTIC.segm, Chromosome != "X", Chromosome != "Y", Chromosome != "M")

write_tsv(GISTIC.segm.no_XYM, "D:/analisi_in_corso/risk_1q_13/GISTIC/GISTIC_1q13_segm_Raw.txt")




#______ create segm_integer file ____________ optional

import <- data.table::fread("D:/analisi_in_corso/risk_1q_13/GISTIC/GISTIC_1q13_segm_Raw.txt")

segm_no_XYM <- filter(segm, chr != "X", chr != "Y", chr != "M")

GISTIC.segm.integer<- data.frame( 
  "Sample" = segm_no_XYM$sample,
  "Chromosome" = segm_no_XYM$chr ,
  "Start Position" = segm_no_XYM$startpos, #+1, # if necessarry manually add 1 base to avoid segment overlapping (GISTIC error)
  "End Position" =segm_no_XYM$endpos ,
  "Num markers" = import$Num.markers ,
  "Seg.CN"= segm_no_XYM$CNvalue_refit_logR
)

write_tsv(GISTIC.segm.integer, "D:/analisi_in_corso/risk_1q_13/GISTIC/GISTIC_1q13_segm_Integer.txt")





