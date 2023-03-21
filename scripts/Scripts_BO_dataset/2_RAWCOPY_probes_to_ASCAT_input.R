# from ASCAT to GISTIC format script

library(data.table)
library(tidyverse)
library(GenomicRanges)
library(foreach)
library(parallel)
library(doParallel)

# ====== data import ======

segmHD1 <- fread("workfiles/ASCAT_refitting/refitted_ASCAT_segments_batchHD1.tsv")
segmHD2 <- fread("workfiles/ASCAT_refitting/refitted_ASCAT_segments_batchHD2.tsv")
segm6.0 <- fread("workfiles/ASCAT_refitting/refitted_ASCAT_segments_batchSNP6.tsv")


#======= Counting of number of probes per segment (required for GISTIC format) ===========

# import PROBES INFO for both types of arrays (pre-generated from rawcopy)
probesHD <- fread("input_data/info_files/probePos_HD.txt")
probes6.0 <- fread("input_data/info_files/probePos_SNP6.txt")


# transform data in Granges
probesHD_GR <- probesHD %>% 
  mutate(start = pos, end = pos) %>% dplyr::rename(chr = chrs ) %>% 
  makeGRangesFromDataFrame(keep.extra.columns = T)

probes6.0_GR <- probes6.0 %>% 
  mutate(start = pos, end = pos) %>% dplyr::rename(chr = chrs ) %>% 
  makeGRangesFromDataFrame(keep.extra.columns = T)

segmHD1_GR <- segmHD1 %>% dplyr::rename(start = startpos, end = endpos) %>% 
  makeGRangesFromDataFrame( keep.extra.columns = T)

segmHD2_GR <- segmHD2 %>% dplyr::rename(start = startpos, end = endpos) %>% 
  makeGRangesFromDataFrame( keep.extra.columns = T)

segm6.0_GR <- segm6.0 %>% dplyr::rename(start = startpos, end = endpos) %>% 
  makeGRangesFromDataFrame( keep.extra.columns = T)


#========== parallel loop to compute number of probes per segment ============

# HD1
registerDoParallel(cores = 8)
result1 <- foreach(i = 1:nrow(segmHD1), .packages = "GenomicRanges", .combine = rbind) %dopar% {
  s <- segmHD1_GR[i,]
  o <- findOverlaps(probesHD_GR, s)
  n <- length(o)
  n
}
segmHD1$num_probes <- result1 %>% as.vector()
stopImplicitCluster()

# HD2
registerDoParallel(cores = 8)
result2 <- foreach(i = 1:nrow(segmHD2), .packages = "GenomicRanges", .combine = rbind) %dopar% {
  s <- segmHD2_GR[i,]
  o <- findOverlaps(probesHD_GR, s)
  n <- length(o)
  n
}
segmHD2$num_probes <- result2 %>% as.vector()
stopImplicitCluster()

# SNP 6.0
registerDoParallel(cores = 8)
result3 <- foreach(i = 1:nrow(segm6.0), .packages = "GenomicRanges", .combine = rbind) %dopar% {
  s <- segm6.0_GR[i,]
  o <- findOverlaps(probes6.0_GR, s)
  n <- length(o)
  n
}
segm6.0$num_probes <- result3 %>% as.vector()
stopImplicitCluster()


#================== merging the 3 batches segments ====================

segm <- rbind(segmHD1, 
              segmHD2, 
              segm6.0)

cor.test(segm$num_probes, segm$length)

segm <- segm %>% mutate(SNP6_type = name %>% str_detect("SNP"))


segm %>% group_by(name) %>% summarise(type = unique(SNP6_type)) %>% View()

segm %>% 
  ggplot(aes(num_probes, length, colour= SNP6_type )) + geom_point(size= 0.5, alpha =0.2) + 
  geom_abline() + geom_abline(intercept = 0, slope = 1350)


#======== integer CN to log2R conversion (required for GISTIC format) ==========

segm$CNvalueRaw_refit_logR <- log2(segm$CNvalueRaw_refit) - 1


segm %>% ggplot(aes(CNvalueRaw_refit_logR)) + geom_density() + xlim(c(-2,2)) + geom_vline(xintercept = 0.585)

segm %>% ggplot(aes(CNvalueRaw_refit)) + geom_density() + xlim(c(0,6)) 


#================ IGV files creation ====================

dir.create("results/IGV_files")

# LogR Raw copy number (refitted for WGD)
IGV_Rawsegm_logR <- data.frame(ID = segm$sample, chrom =segm$chr, loc.start= segm$startpos,
                                loc.end = segm$endpos, num.mark	= segm$num_probes,
                                seg.mean = segm$CNvalueRaw_refit_logR)

readr::write_tsv(IGV_Rawsegm_logR, "results/IGV_files/IGV_BO_Rawsegm_logR.seg" )


# Integer Raw copy number (refitted for WGD)
IGV_Rawsegm <- data.frame(ID = segm$sample, chrom =segm$chr, loc.start= segm$startpos,
                               loc.end = segm$endpos, num.mark	= segm$num_probes,
                               seg.mean = segm$CNvalueRaw_refit)

readr::write_tsv(IGV_Rawsegm, "results/IGV_files/IGV_BO_Rawsegm.seg" )



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

write_tsv(GISTIC.segm.no_XYM, "workfiles/GISTIC_1q13_segm_Raw.txt")



