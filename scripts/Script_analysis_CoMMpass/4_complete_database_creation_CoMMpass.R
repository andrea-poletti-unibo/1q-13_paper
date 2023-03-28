
library(tidyverse)
library(data.table)

#_____import the txt file from portal with CoMMpass clinical and FISH data______

FISH_CLIN <- fread("input_data/clinical_data_COMMPASS-IA13_210122.txt")

#_____import the prepared CNA calls of DELs and AMPs from GISTIC run (broad and focal)______

#broad maj
AMPcalls_broad_maj <- fread("workfiles/BROAD_CALLS/CoMMpass/AMP_CALLS_major_broad.txt")
DELcalls_broad_maj <- fread("workfiles/BROAD_CALLS/CoMMpass/DEL_CALLS_major_broad.txt")
setnames(AMPcalls_broad_maj, "V1", "sample")
setnames(DELcalls_broad_maj, "V1", "sample")


# focal maj
AMPcalls_focal_maj <- fread("workfiles/FOCAL_CALLS/CoMMpass/AMP_CALLS_major_focal_Compass.txt")
DELcalls_focal_maj <- fread("workfiles/FOCAL_CALLS/CoMMpass//DEL_CALLS_major_focal_Compass.txt")
setnames(AMPcalls_focal_maj, "V1", "sample")
setnames(DELcalls_focal_maj, "V1", "sample")


#==================== creation of complete database =========================

complete_callset <- full_join(AMPcalls_broad_maj, 
                              DELcalls_broad_maj, by = "sample") %>%
  full_join(AMPcalls_focal_maj, by = "sample") %>% 
  full_join(DELcalls_focal_maj, by = "sample")

complete_database <- inner_join(FISH_CLIN, complete_callset, by = c("Study_Visit_iD"="sample"))


#_____ hyperdiploidy calculation ________

# ALL level of clonality

complete_database$HD_chr_3  <- ifelse( with(complete_database, `AMP_maj_broad_chr_3p`==1 | `AMP_maj_broad_chr_3q`==1 ), 1,0) 
complete_database$HD_chr_5  <- ifelse( with(complete_database, `AMP_maj_broad_chr_5p`==1 | `AMP_maj_broad_chr_5q`==1 ), 1,0) 
complete_database$HD_chr_7  <- ifelse( with(complete_database, `AMP_maj_broad_chr_7p`==1 | `AMP_maj_broad_chr_7q`==1 ), 1,0) 
complete_database$HD_chr_9  <- ifelse( with(complete_database, `AMP_maj_broad_chr_9p`==1 | `AMP_maj_broad_chr_9q`==1 ), 1,0) 
complete_database$HD_chr_11 <- ifelse( with(complete_database, `AMP_maj_broad_chr_11p`==1| `AMP_maj_broad_chr_11q`==1 ), 1,0) 
complete_database$HD_chr_15 <- complete_database$`AMP_maj_broad_chr_15q`
complete_database$HD_chr_19 <- ifelse( with(complete_database, `AMP_maj_broad_chr_19p`==1 | `AMP_maj_broad_chr_19q`==1 ), 1,0) 
complete_database$HD_chr_21 <- complete_database$`AMP_maj_broad_chr_21q`

complete_database$HD.count <- with(complete_database, HD_chr_3+HD_chr_5+HD_chr_7+HD_chr_9+HD_chr_11+HD_chr_15+HD_chr_19+HD_chr_21) 

complete_database$HD.count %>% table

complete_database$HyperDiploidy <- ifelse(complete_database$HD.count>=2,1,0) 

complete_database$HyperDiploidy %>% table




#______________ Compute Traslocation_IgH variable ________________________

complete_database$T_IgH <- ifelse( complete_database$SeqWGS_WHSC1_CALL + 
                                     complete_database$SeqWGS_CCND1_CALL + 
                                     complete_database$SeqWGS_CCND3_CALL + 
                                     complete_database$SeqWGS_MAF_CALL + 
                                     complete_database$SeqWGS_MAFB_CALL > 0, 1, 0)

table(complete_database$T_IgH)

table(complete_database$T_IgH, complete_database$HyperDiploidy)


complete_database$t_others <- complete_database$SeqWGS_CCND3_CALL + complete_database$SeqWGS_MAF_CALL + complete_database$SeqWGS_MAFB_CALL


#____ Build CALLS for arms with broad+focal_____

#amps
complete_database$AMP_maj_call_chr_1q <- ifelse( complete_database$AMP_maj_broad_chr_1q == 1 | complete_database$AMP_maj_focal_ANP32E== 1 | complete_database$AMP_maj_focal_CKS1B ==1 | complete_database$AMP_maj_focal_MCL1 ==1, 1, 0)
complete_database$AMP_maj_call_chr_8q <- ifelse( complete_database$AMP_maj_broad_chr_8q == 1 | complete_database$AMP_maj_focal_MYC == 1, 1, 0)

#dels
complete_database$DEL_maj_call_chr_17p <- ifelse( complete_database$DEL_maj_broad_chr_17p == 1 | complete_database$DEL_maj_focal_TP53 == 1, 1, 0)
complete_database$DEL_maj_call_chr_1p <-  ifelse( complete_database$DEL_maj_broad_chr_1p == 1 |  complete_database$DEL_maj_focal_CDKN2C == 1 | complete_database$DEL_maj_focal_FAM46C == 1 | complete_database$DEL_maj_focal_FAF1 == 1,  1, 0)
complete_database$DEL_maj_call_chr_13q <- ifelse( complete_database$DEL_maj_broad_chr_13q == 1 | complete_database$DEL_maj_focal_RB1 == 1, 1, 0)
complete_database$DEL_maj_call_chr_14q <- ifelse( complete_database$DEL_maj_broad_chr_14q == 1 | complete_database$DEL_maj_focal_TRAF3 == 1, 1, 0)
complete_database$DEL_maj_call_chr_16q <- ifelse( complete_database$DEL_maj_broad_chr_16q == 1 | complete_database$DEL_maj_focal_CYLD == 1, 1, 0)


#___________ Build 1q&13 classifier labels _____________
complete_database$AMP_1q_genes_all <- with(complete_database, 
                                ifelse(`AMP_maj_focal_ANP32E`==1 | `AMP_maj_focal_MCL1` ==1 | `AMP_maj_focal_CKS1B`==1, 1,0 ))

complete_database$MMrisk_1q_all <- ifelse( complete_database$`AMP_maj_broad_chr_1q` == 1 | complete_database$AMP_1q_genes_all ==1, 1, 0)
complete_database$MMrisk_13_all <- ifelse( complete_database$`DEL_maj_broad_chr_13q` ==1 | complete_database$`DEL_maj_focal_RB1` ==1, 1, 0)

complete_database$MMrisk_CLASSIFIER_ALL <- 3  - complete_database$MMrisk_1q_all - complete_database$MMrisk_13_all
complete_database$MMrisk_CLASSIFIER_ALL <- dplyr::recode(as.character(complete_database$MMrisk_CLASSIFIER_ALL), "1"="risk1","2"="risk2","3"="risk3")

complete_database$MMrisk_CLASSIFIER_ALL %>% table


#_____ export ______
write_tsv(complete_database, "D:/analisi_in_corso/risk_1q_13/CoMMpass/complete_database_1q_13_CoMMpass_220121.txt")

write_tsv(complete_database, "c:/Users/mm_gr/Alma Mater Studiorum Università di Bologna/PROJECT 1q & 13 - Documenti/complete_database_1q_13_CoMMpass_220121.txt")





