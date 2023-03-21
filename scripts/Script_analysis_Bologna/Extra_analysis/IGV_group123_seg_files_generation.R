
library(data.table)
library(tidyverse)
library(RODBC)



db_clin <- RODBC::odbcConnectAccess2007("C:/Users/andre/Alma Mater Studiorum Università di Bologna/MM group - Documenti/EVERGREEN_inter-projects_files/MM_db_clinical_records.accdb")
clin <- sqlFetch(db_clin, "Extraction_1q13_301120", dec=",", na.strings="nv" )

db_exp<- RODBC::odbcConnectAccess2007("C:/Users/andre/Alma Mater Studiorum Università di Bologna/MM group - Documenti/EVERGREEN_inter-projects_files/MM_db_experiments.accdb")
calls <- sqlFetch(db_exp, "SNP_calls_1q_13_181019", dec=",")

import <- inner_join(clin, calls, by=c("SNP_ARRAY_ESORDIO"="SNP"))

import$PROTOCOL_REV <- recode(import$PROTOCOL, "FP-EMN02"="AMB", "FP - BO2005"="AMB", "DARA"="AMB", "FORTE"="AMB")
import$IG_ISOTYPE_REV <- recode(import$IG_ISOTYPE, "IgA+BJ"="IgA", "IgA/BJ"="IgA")


#______ rename variables _______
names(import)
import <- import %>% rename(ASCT = TX_I_LINE_1_0,
                            MAINTENANCE = MAINTENANCE_YES_NO,
                            CONSOLIDATION = CONSOLIDATION_YES_NO)


#______ MM RISK CREATION _____
import$AMP_1q_genes_all <- with(import, 
                                ifelse(`AMP_maj-focal_ANP32E`==1 | `AMP_maj-focal_MCL1` ==1 | `AMP_maj-focal_CKS1B`==1, 1,0 ))

import$MMrisk_1q_all <- ifelse( import$`AMP_maj-broad_chr_1q` == 1 | import$AMP_1q_genes_all ==1, 1, 0)
import$MMrisk_13_all <- ifelse( import$`DEL_maj-broad_chr_13q` ==1 | import$`DEL_maj-focal_RB1` ==1, 1, 0)

import$MMrisk_CLASS <- 3  - import$MMrisk_1q_all - import$MMrisk_13_all
import$MMrisk_class <- paste0("risk_",import$MMrisk_CLASS) %>% as.factor() %>% relevel("risk_2")

import$MMrisk1 <- ifelse(import$MMrisk_CLASS==1, 1, 0)
import$MMrisk2 <- ifelse(import$MMrisk_CLASS==2, 1, 0)
import$MMrisk3 <- ifelse(import$MMrisk_CLASS==3, 1, 0)


SNP_RISK <- import %>% select(SNP_FILE_NAME, MMrisk_class)

IGV_segs <- fread("IGV/IGV_Rawsegm_medianCentered_Scaled.seg")
IGV_segs_2 <- IGV_segs

IGV_segs_2$ID <- IGV_segs$ID %>% str_remove("^X") %>% 
  str_replace("\\.GenomeWideSNP_6\\.", "\\(GenomeWideSNP_6\\)") %>% 
  str_replace("\\.CytoScanHD_Array\\.", "\\(CytoScanHD_Array\\)")

IGV_segs_2$ID %>% unique()


SNP_RISK$SNP_FILE_NAME %in% IGV_segs_2$ID %>% table

R1_IGV <- IGV_segs_2[IGV_segs_2$ID %in% SNP_RISK$SNP_FILE_NAME[SNP_RISK$MMrisk_class=="risk_1"] , ]
R1_IGV$ID %>% unique() %>% length()


R2_IGV <- IGV_segs_2[IGV_segs_2$ID %in% SNP_RISK$SNP_FILE_NAME[SNP_RISK$MMrisk_class=="risk_2"] , ]
R2_IGV$ID %>% unique() %>% length()


R3_IGV <- IGV_segs_2[IGV_segs_2$ID %in% SNP_RISK$SNP_FILE_NAME[SNP_RISK$MMrisk_class=="risk_3"] , ]
R3_IGV$ID %>% unique() %>% length()

131 + 170 + 211

write_tsv(R1_IGV, "IGV/IGV_Rawsegm_medianCentered_Scaled_MMrisk1.seg")
write_tsv(R2_IGV, "IGV/IGV_Rawsegm_medianCentered_Scaled_MMrisk2.seg")
write_tsv(R3_IGV, "IGV/IGV_Rawsegm_medianCentered_Scaled_MMrisk3.seg")

