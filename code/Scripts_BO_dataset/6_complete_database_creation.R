# script to create the complete database (CNA clonal broad + FISH + clinical)

library(tidyverse)
library(data.table)

#_____import the data with FISH and clincal data______

FISH_CLIN <- fread("data/clinical_data_BOLO_210122_anonim.txt")


#_____import the prepared CNA calls of DELs and AMPs from GISTIC run (broad and focal)______

#broad maj
AMPcalls_broad_maj <- fread("workfiles/BROAD_CALLS/AMP_CALLS_major_broad.txt")
DELcalls_broad_maj <- fread("workfiles/BROAD_CALLS/DEL_CALLS_major_broad.txt")
setnames(AMPcalls_broad_maj, "V1", "sample")
setnames(DELcalls_broad_maj, "V1", "sample")

#broad all
AMPcalls_broad_all <- fread("workfiles/BROAD_CALLS/AMP_CALLS_all_broad.txt")
DELcalls_broad_all <- fread("workfiles/BROAD_CALLS/DEL_CALLS_all_broad.txt")
setnames(AMPcalls_broad_all, "V1", "sample")
setnames(DELcalls_broad_all, "V1", "sample")



# focal maj
AMPcalls_focal_maj <- fread("workfiles/FOCAL_CALLS/AMP_CALLS_major_focal.txt")
DELcalls_focal_maj <- fread("workfiles/FOCAL_CALLS/DEL_CALLS_major_focal.txt")
setnames(AMPcalls_focal_maj, "V1", "sample")
setnames(DELcalls_focal_maj, "V1", "sample")

# focal all
AMPcalls_focal_all <- fread("workfiles/FOCAL_CALLS/AMP_CALLS_all_focal.txt")
DELcalls_focal_all <- fread("workfiles/FOCAL_CALLS/DEL_CALLS_all_focal.txt")
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

complete_database <- left_join(FISH_CLIN, complete_callset, by = c("SNP_FILE_NAME"="sample"))


table(complete_callset$sample %in% FISH_CLIN$SNP_FILE_NAME)

complete_callset$sample[!complete_callset$sample %in% FISH_CLIN$SNP_FILE_NAME]



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

#_________ trasnslocations calls ________

complete_database <- complete_database %>% dplyr::mutate(t_4_14 = FISH_T_4_14,
                                                         t_11_14 = FISH_T_11_14,
                                                         t_6_14 = FISH_T_6_14,
                                                         t_14_16 = FISH_T_14_16,
                                                         t_14_20 = FISH_T_14_20)

#______ Build CALLS for arms with broad+focal_____

#amps
complete_database$AMP_maj_call_chr_1q <- ifelse( complete_database$AMP_maj_broad_chr_1q == 1 | complete_database$AMP_maj_focal_ANP32E== 1 | complete_database$AMP_maj_focal_CKS1B ==1 | complete_database$AMP_maj_focal_MCL1 ==1, 1, 0)
complete_database$AMP_maj_call_chr_8q <- ifelse( complete_database$AMP_maj_broad_chr_8q == 1 | complete_database$AMP_maj_focal_MYC == 1, 1, 0)
#dels
complete_database$DEL_maj_call_chr_17p <- ifelse( complete_database$DEL_maj_broad_chr_17p == 1 | complete_database$DEL_maj_focal_TP53 == 1, 1, 0)
complete_database$DEL_maj_call_chr_1p <- ifelse( complete_database$DEL_maj_broad_chr_1p == 1 | complete_database$DEL_maj_focal_CDKN2C == 1 | complete_database$DEL_maj_focal_FAM46C == 1 | complete_database$DEL_maj_focal_FAF1 == 1,  1, 0)
complete_database$DEL_maj_call_chr_13q <- ifelse( complete_database$DEL_maj_broad_chr_13q == 1 | complete_database$DEL_maj_focal_RB1 == 1, 1, 0)
complete_database$DEL_maj_call_chr_14q <- ifelse( complete_database$DEL_maj_broad_chr_14q == 1 | complete_database$DEL_maj_focal_TRAF3 == 1, 1, 0)
complete_database$DEL_maj_call_chr_16q <- ifelse( complete_database$DEL_maj_broad_chr_16q == 1 | complete_database$DEL_maj_focal_CYLD == 1, 1, 0)


#_____ export ______
write_tsv(complete_database, "results/complete_database_BO_1q_13.txt")


