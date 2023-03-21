# script to create the complete database (CNA clonal broad + FISH + clinical)

library(tidyverse)
library(data.table)

setwd("D:/analisi_in_corso/risk_1q_13/")

#_____import the excel (from Cinzia) with FISH and clincal data______

FISH_CLIN <- readxl::read_xlsx("D:/analisi_in_corso/risk_1q_13/1q_13_clinicaldata_310519.xlsx")

#_____import the prepared CNA calls of DELs and AMPs from GISTIC run (broad and focal)______

#broad maj
AMPcalls_broad_maj <- fread("GISTIC/BROAD_CALLS/AMP_CALLS_major_broad.txt")
DELcalls_broad_maj <- fread("GISTIC/BROAD_CALLS/DEL_CALLS_major_broad.txt")
setnames(AMPcalls_broad_maj, "V1", "sample")
setnames(DELcalls_broad_maj, "V1", "sample")

#broad all
AMPcalls_broad_all <- fread("GISTIC/BROAD_CALLS/AMP_CALLS_all_broad.txt")
DELcalls_broad_all <- fread("GISTIC/BROAD_CALLS/DEL_CALLS_all_broad.txt")
setnames(AMPcalls_broad_all, "V1", "sample")
setnames(DELcalls_broad_all, "V1", "sample")



# focal maj
AMPcalls_focal_maj <- fread("GISTIC/FOCAL_CALLS/AMP_CALLS_major_focal.txt")
DELcalls_focal_maj <- fread("GISTIC/FOCAL_CALLS/DEL_CALLS_major_focal.txt")
setnames(AMPcalls_focal_maj, "V1", "sample")
setnames(DELcalls_focal_maj, "V1", "sample")

# focal all
AMPcalls_focal_all <- fread("GISTIC/FOCAL_CALLS/AMP_CALLS_all_focal.txt")
DELcalls_focal_all <- fread("GISTIC/FOCAL_CALLS/DEL_CALLS_all_focal.txt")
setnames(AMPcalls_focal_all, "V1", "sample")
setnames(DELcalls_focal_all, "V1", "sample")



# correction of samplenames
correctnames<- AMPcalls_broad_maj$sample %>% str_replace("^X","") %>% str_replace("\\.CytoScanHD_Array\\.","(CytoScanHD_Array)") %>% str_replace("\\.GenomeWideSNP_6\\.","(GenomeWideSNP_6)")  %>% str_replace("\\.1035_"," 1035_")   

correctnames %in% FISH_CLIN$nome_CEL_sample %>% table
correctnames[!(correctnames %in% FISH_CLIN$nome_CEL_sample)] #sample MM_cytoscanHD_08282017_820_(CytoScanHD_Array) now doesnt'have clinical data! (was in "1q_13_clinicaldata_240519.xlsx" but not in "1q_13_clinicaldata_310519.xlsx" ) 


AMPcalls_broad_maj$sample <- AMPcalls_broad_maj$sample %>% str_replace("^X","") %>% str_replace("\\.CytoScanHD_Array\\.","(CytoScanHD_Array)") %>% str_replace("\\.GenomeWideSNP_6\\.","(GenomeWideSNP_6)")  %>% str_replace("\\.1035_"," 1035_")   
DELcalls_broad_maj$sample <- DELcalls_broad_maj$sample %>% str_replace("^X","") %>% str_replace("\\.CytoScanHD_Array\\.","(CytoScanHD_Array)") %>% str_replace("\\.GenomeWideSNP_6\\.","(GenomeWideSNP_6)")  %>% str_replace("\\.1035_"," 1035_")   

AMPcalls_broad_all$sample <- AMPcalls_broad_all$sample %>% str_replace("^X","") %>% str_replace("\\.CytoScanHD_Array\\.","(CytoScanHD_Array)") %>% str_replace("\\.GenomeWideSNP_6\\.","(GenomeWideSNP_6)")  %>% str_replace("\\.1035_"," 1035_")   
DELcalls_broad_all$sample <- DELcalls_broad_all$sample %>% str_replace("^X","") %>% str_replace("\\.CytoScanHD_Array\\.","(CytoScanHD_Array)") %>% str_replace("\\.GenomeWideSNP_6\\.","(GenomeWideSNP_6)")  %>% str_replace("\\.1035_"," 1035_")   


AMPcalls_focal_maj$sample <- AMPcalls_focal_maj$sample %>% str_replace("^X","") %>% str_replace("\\.CytoScanHD_Array\\.","(CytoScanHD_Array)") %>% str_replace("\\.GenomeWideSNP_6\\.","(GenomeWideSNP_6)")  %>% str_replace("\\.1035_"," 1035_")   
DELcalls_focal_maj$sample <- DELcalls_focal_maj$sample %>% str_replace("^X","") %>% str_replace("\\.CytoScanHD_Array\\.","(CytoScanHD_Array)") %>% str_replace("\\.GenomeWideSNP_6\\.","(GenomeWideSNP_6)")  %>% str_replace("\\.1035_"," 1035_")   

AMPcalls_focal_all$sample <- AMPcalls_focal_all$sample %>% str_replace("^X","") %>% str_replace("\\.CytoScanHD_Array\\.","(CytoScanHD_Array)") %>% str_replace("\\.GenomeWideSNP_6\\.","(GenomeWideSNP_6)")  %>% str_replace("\\.1035_"," 1035_")   
DELcalls_focal_all$sample <- DELcalls_focal_all$sample %>% str_replace("^X","") %>% str_replace("\\.CytoScanHD_Array\\.","(CytoScanHD_Array)") %>% str_replace("\\.GenomeWideSNP_6\\.","(GenomeWideSNP_6)")  %>% str_replace("\\.1035_"," 1035_")   


#==================== creation of complete database =========================

complete_callset <- full_join(AMPcalls_broad_maj, DELcalls_broad_maj, by = "sample") %>% 
  full_join(AMPcalls_broad_all, by = "sample") %>% 
  full_join(DELcalls_broad_all, by = "sample") %>% 
  full_join(AMPcalls_focal_maj, by = "sample") %>% 
  full_join(DELcalls_focal_maj, by = "sample") %>% 
  full_join(AMPcalls_focal_all, by = "sample") %>% 
  full_join(DELcalls_focal_all, by = "sample") 
  
complete_database <- left_join(FISH_CLIN, complete_callset, by = c("nome_CEL_sample"="sample"))


#_____ hyperdiploidy calculation ________

# ALL level of clonality

complete_database$HD_chr_3  <- ifelse( with(complete_database, `AMP_all_broad_chr_3p`==1 | `AMP_all_broad_chr_3q`==1 ), 1,0) 
complete_database$HD_chr_5  <- ifelse( with(complete_database, `AMP_all_broad_chr_5p`==1 | `AMP_all_broad_chr_5q`==1 ), 1,0) 
complete_database$HD_chr_7  <- ifelse( with(complete_database, `AMP_all_broad_chr_7p`==1 | `AMP_all_broad_chr_7q`==1 ), 1,0) 
complete_database$HD_chr_9  <- ifelse( with(complete_database, `AMP_all_broad_chr_9p`==1 | `AMP_all_broad_chr_9q`==1 ), 1,0) 
complete_database$HD_chr_11 <- ifelse( with(complete_database, `AMP_all_broad_chr_11p`==1| `AMP_all_broad_chr_11q`==1 ), 1,0) 
complete_database$HD_chr_15 <- complete_database$`AMP_all_broad_chr_15q`
complete_database$HD_chr_19 <- ifelse( with(complete_database, `AMP_all_broad_chr_19p`==1 | `AMP_all_broad_chr_19q`==1 ), 1,0) 
complete_database$HD_chr_21 <- complete_database$`AMP_all_broad_chr_21q`

complete_database$HD.count <- with(complete_database, HD_chr_3+HD_chr_5+HD_chr_7+HD_chr_9+HD_chr_11+HD_chr_15+HD_chr_19+HD_chr_21) 

complete_database$HD.count %>% table

complete_database$HyperDiploidy <- ifelse(complete_database$HD.count>=2,1,0) 

complete_database$HyperDiploidy %>% table

#_____ export ______
write_tsv(complete_database, "complete_database_1q_13_181019.txt")


