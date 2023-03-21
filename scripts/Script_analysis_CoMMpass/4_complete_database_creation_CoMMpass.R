# script to create the complete database (CNA clonal broad + FISH + clinical)

library(tidyverse)
library(data.table)

setwd("D:/analisi_in_corso/risk_1q_13/")

#_____import the excel (from Cinzia) with FISH and clincal data______

FISH_CLIN <- fread("D:/analisi_in_corso/risk_1q_13/CoMMpass/Traslocations_and_survival.txt")

#_____import the prepared CNA calls of DELs and AMPs from GISTIC run (broad and focal)______

#broad maj
AMPcalls_broad_maj <- fread("CoMMpass/BROAD_FOCAL_CALLS/CoMMpass_AMP_maj_broad.txt")
DELcalls_broad_maj <- fread("CoMMpass/BROAD_FOCAL_CALLS/CoMMpass_DEL_maj_broad.txt")
setnames(AMPcalls_broad_maj, "V1", "sample")
setnames(DELcalls_broad_maj, "V1", "sample")

#broad all
AMPcalls_broad_all <- fread("CoMMpass/BROAD_FOCAL_CALLS/CoMMpass_AMP_all_broad.txt")
DELcalls_broad_all <- fread("CoMMpass/BROAD_FOCAL_CALLS/CoMMpass_DEL_all_broad.txt")
setnames(AMPcalls_broad_all, "V1", "sample")
setnames(DELcalls_broad_all, "V1", "sample")


# focal maj
AMPcalls_focal_maj <- fread("CoMMpass/BROAD_FOCAL_CALLS/AMP_CALLS_major_focal_Compass.txt")
DELcalls_focal_maj <- fread("CoMMpass/BROAD_FOCAL_CALLS/DEL_CALLS_major_focal_Compass.txt")
setnames(AMPcalls_focal_maj, "V1", "sample")
setnames(DELcalls_focal_maj, "V1", "sample")

# focal all
AMPcalls_focal_all <- fread("CoMMpass/BROAD_FOCAL_CALLS/AMP_CALLS_all_focal.txt")
DELcalls_focal_all <- fread("CoMMpass/BROAD_FOCAL_CALLS/DEL_CALLS_all_focal.txt")
setnames(AMPcalls_focal_all, "V1", "sample")
setnames(DELcalls_focal_all, "V1", "sample")



#==================== creation of complete database =========================

complete_callset <- full_join(AMPcalls_broad_maj, DELcalls_broad_maj, by = "sample") %>% 
  full_join(AMPcalls_broad_all, by = "sample") %>% 
  full_join(DELcalls_broad_all, by = "sample") %>%   
  full_join(AMPcalls_focal_maj, by = "sample") %>% 
  full_join(DELcalls_focal_maj, by = "sample") %>% 
  full_join(AMPcalls_focal_all, by = "sample") %>% 
  full_join(DELcalls_focal_all, by = "sample") 

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

#_____ export ______
write_tsv(complete_database, "D:/analisi_in_corso/risk_1q_13/CoMMpass/complete_database_1q_13_CoMMpass_220121.txt")

write_tsv(complete_database, "c:/Users/mm_gr/Alma Mater Studiorum Università di Bologna/PROJECT 1q & 13 - Documenti/complete_database_1q_13_CoMMpass_220121.txt")



################ Complete_pt2 : add clinical infos #################

library(tidyverse)
library(data.table)

import <- fread("c:/Users/andre/Alma Mater Studiorum Università di Bologna/PROJECT 1q & 13 - Documenti/complete_database_1q_13_CoMMpass_220121.txt")

import$PUBLIC_ID <- import$Study_Visit_iD %>% str_extract("MMRF_[0-9]+")


pervisit <-fread("C:/Users/andre/Alma Mater Studiorum Università di Bologna/MM group - Documenti/EVERGREEN_inter-projects_files/COMMPASS/IA13a/Clinical_Flat_Files/CoMMpass_IA13_FlatFiles/MMRF_CoMMpass_IA13_PER_PATIENT_VISIT.csv")

baseline <- pervisit %>% filter(VJ_INTERVAL=="Baseline", grepl("_1$",pervisit$SPECTRUM_SEQ))
duplicated(baseline$PUBLIC_ID) %>% sum

res <- left_join(import, baseline, by="PUBLIC_ID")

#_____ export ______
write_tsv(res, "D:/analisi_in_corso/risk_1q_13/CoMMpass/complete_database_1q_13_CoMMpass_240121.txt")

write_tsv(res, "C:/Users/andre/Alma Mater Studiorum Università di Bologna/PROJECT 1q & 13 - Documenti/complete_database_1q_13_CoMMpass_240121.txt")



