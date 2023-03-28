
####### from Commpass IA13a WGS segmensts (portal download) to GISTIC format ###########

#__________ load libraries __________

library(data.table)
library(tidyverse)

#__________ data import e prep__________

import <- fread("input_data/CoMMpass/MMRF_CoMMpass_IA13a_CNA_LongInsert.seg") # load Compass Segment data 

# filter per first visit (1_BM in sample name)

segm <- import %>% filter(ID %>% str_detect("1_BM"))

segm$ID %>% unique
segm$ID %>% unique %>% length() # 870 samples baseline BM

#__________ creation of a dataframe in GISTIC format___________

# Change the name of columns

GISTIC.segm<- data.frame( 
                       "Sample" = segm$ID , 
                       "Chromosome" = segm$chrom , 
                       "Start Position" = segm$loc.start , 
                       "End Position" =segm$loc.end , 
                       "Num markers" = segm$num.mark , 
                       "Seg.CN"= segm$seg.mean 
                       )

# filter to exclude chr X, Y and Mitocondrial
GISTIC.segm.no_XYM<- filter(GISTIC.segm, Chromosome != "X", Chromosome != "Y", Chromosome != "M")


# save output
write_tsv(GISTIC.segm.no_XYM, "workfiles/GISTIC_1q13_segm_CoMMpass.txt") 





#================ IGV files creation ====================

dir.create("results/IGV_files/CoMMpass", recursive = T, showWarnings = F)

segm$CN <- 2^(segm$seg.mean + 1)

segm %>% ggplot(aes(CN)) + geom_density() + xlim(0,4.5)

# LogR Raw copy number (refitted for WGD)
IGV_Rawsegm_logR <- data.frame(ID = segm$ID, chrom =segm$chrom, loc.start= segm$loc.start,
                               loc.end = segm$loc.end, num.mark	= segm$num.mark,
                               seg.mean = segm$seg.mean)

readr::write_tsv(IGV_Rawsegm_logR, "results/IGV_files/CoMMpass/IGV_CoMM_logR.seg" )


# Integer Raw copy number (refitted for WGD)
IGV_Rawsegm <- data.frame(ID = segm$ID, chrom =segm$chrom, loc.start= segm$loc.start,
                          loc.end = segm$loc.end, num.mark	= segm$num.mark,
                          seg.mean = segm$CN)

readr::write_tsv(IGV_Rawsegm, "results/IGV_files/CoMMpass/IGV_CoMM.seg" )







