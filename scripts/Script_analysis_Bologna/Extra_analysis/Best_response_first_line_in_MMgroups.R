library(data.table)
library(tidyverse)
library(RODBC)


db_clin <- RODBC::odbcConnectAccess2007("C:/Users/mm_gr/Alma Mater Studiorum Università di Bologna/MM group - Documenti/EVERGREEN_inter-projects_files/MM_db_clinical_records.accdb")
clin <- sqlFetch(db_clin, "Extraction_1q13_301120", dec=",", na.strings="nv" )

db_exp<- RODBC::odbcConnectAccess2007("C:/Users/mm_gr/Alma Mater Studiorum Università di Bologna/MM group - Documenti/EVERGREEN_inter-projects_files/MM_db_experiments.accdb")
calls <- sqlFetch(db_exp, "SNP_calls_1q_13_181019", dec=",")

import <- inner_join(clin, calls, by=c("SNP_ARRAY_ESORDIO"="SNP"))

import$IG_ISOTYPE_REV <- recode(import$IG_ISOTYPE, "IgA+BJ"="IgA", "IgA/BJ"="IgA")

import$PROTOCOL_REV <- recode(import$PROTOCOL, "FP-EMN02"="AMB", "FP - BO2005"="AMB", "DARA"="AMB", "FORTE"="AMB")


#______ MM RISK CREATION _____
import$AMP_1q_genes_all <- with(import, 
                                ifelse(`AMP_maj-focal_ANP32E`==1 | `AMP_maj-focal_MCL1` ==1 | `AMP_maj-focal_CKS1B`==1, 1,0 ))

import$MMrisk_1q_all <- ifelse( import$`AMP_maj-broad_chr_1q` == 1 | import$AMP_1q_genes_all ==1, 1, 0)
import$MMrisk_13_all <- ifelse( import$`DEL_maj-broad_chr_13q` ==1 | import$`DEL_maj-focal_RB1` ==1, 1, 0)

import$MMrisk_CLASSIFIER_ALL <- 3  - import$MMrisk_1q_all - import$MMrisk_13_all
import$MMrisk_class <- paste0("risk_",import$MMrisk_CLASSIFIER_ALL) %>% as.factor() %>% relevel("risk_2")

import$MMrisk_AB_ALL <- ifelse(import$MMrisk_CLASSIFIER_ALL == 2 & import$MMrisk_1q_all ==1,
                               "2_1q", 
                               ifelse(import$MMrisk_CLASSIFIER_ALL == 2 & import$MMrisk_13_all ==1,
                                      "2_13",
                                      import$MMrisk_CLASSIFIER_ALL))

import$risk1_VS_1q_VS_other <- recode(import$MMrisk_AB_ALL, "1"="1q&13", "2_1q"="1q_only", "2_13"="other", "3"="other")


import$MMrisk1 <- ifelse(import$MMrisk_CLASSIFIER_ALL==1, 1, 0)
import$MMrisk2 <- ifelse(import$MMrisk_CLASSIFIER_ALL==2, 1, 0)
import$MMrisk3 <- ifelse(import$MMrisk_CLASSIFIER_ALL==3, 1, 0)



import$ULTRA_MMrisk <- ifelse(import$ISS == 3 & import$MMrisk_CLASSIFIER_ALL == 1,
                              "Ultra_High",
                              ifelse(import$ISS==1 & import$MMrisk_CLASSIFIER_ALL == 3,
                                     "Ultra_Low",
                                     "Other"))

import$ULTRA_High <- ifelse(import$ULTRA_MMrisk=="Ultra_High",1,0)
import$ULTRA_Low  <- ifelse(import$ULTRA_MMrisk=="Ultra_Low",1,0)



#_________________ confront best clinical response ________________

import$FIRST_LINE_BEST_RESPONSE %>% table

import$FL_BR_more_nCR <- ifelse(import$FIRST_LINE_BEST_RESPONSE %in% c("nCR", "CR", "sCR"), 1,0)
import$FL_BR_more_nCR %>% table

import$FL_BR_more_VGPR <- ifelse(import$FIRST_LINE_BEST_RESPONSE %in% c("nCR", "CR", "sCR", "VGPR"), 1,0)
import$FL_BR_more_VGPR %>% table


gmodels::CrossTable(import$MMrisk_class, import$FL_BR_more_nCR, fisher = T, prop.chisq = F, prop.t = F)

gmodels::CrossTable(import$MMrisk_class, import$FL_BR_more_VGPR, fisher = T, prop.chisq = F, prop.t = F)








