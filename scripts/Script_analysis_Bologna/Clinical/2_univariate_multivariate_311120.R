library(data.table)
library(tidyverse)
library(survival)
library(survminer)
library(RODBC)

db_clin <- RODBC::odbcConnectAccess2007("C:/Users/mm_gr/Alma Mater Studiorum Università di Bologna/MM group - Documenti/EVERGREEN_inter-projects_files/MM_db_clinical_records.accdb")
sqlTables(db_clin)
clin <- sqlFetch(db_clin, "Extraction_1q13_301120", dec=",", na.strings="nv" )


db_exp<- RODBC::odbcConnectAccess2007("C:/Users/mm_gr/Alma Mater Studiorum Università di Bologna/MM group - Documenti/EVERGREEN_inter-projects_files/MM_db_experiments.accdb")
sqlTables(db_exp)
calls <- sqlFetch(db_exp, "SNP_calls_1q_13_181019", dec=",")


calls$SNP
clin$SNP_ARRAY_ESORDIO

import <- inner_join(clin, calls, by=c("SNP_ARRAY_ESORDIO"="SNP"))


calls$SNP[!(calls$SNP %in% clin$SNP_ARRAY_ESORDIO)]
clin$SNP_ARRAY_ESORDIO[!(clin$SNP_ARRAY_ESORDIO %in% calls$SNP)]



import$PROTOCOL %>% table

import$PROTOCOL_REV <- recode(import$PROTOCOL, "FP-EMN02"="AMB", "FP - BO2005"="AMB", "DARA"="AMB", "FORTE"="AMB")

import$PROTOCOL_REV %>% table

import$IG_ISOTYPE_REV <- recode(import$IG_ISOTYPE, "IgA+BJ"="IgA", "IgA/BJ"="IgA")

import$IG_ISOTYPE_REV  %>% table


import$NUMBERS_TX_FRONT_LINE %>% table

#______ MM RISK CREATION _____
import$AMP_1q_genes_all <- with(import, 
                                ifelse(`AMP_maj-focal_ANP32E`==1 | `AMP_maj-focal_MCL1` ==1 | `AMP_maj-focal_CKS1B`==1, 1,0 ))


table(import$`AMP_maj-broad_chr_1q`, import$AMP_1q_genes_all)
table(import$`DEL_maj-broad_chr_13q`, import$`DEL_maj-focal_RB1`)


import$MMrisk_1q_all <- ifelse( import$`AMP_maj-broad_chr_1q` == 1 | import$AMP_1q_genes_all ==1, 1, 0)
import$MMrisk_13_all <- ifelse( import$`DEL_maj-broad_chr_13q` ==1 | import$`DEL_maj-focal_RB1` ==1, 1, 0)

import$MMrisk_CLASSIFIER_ALL <- 3  - import$MMrisk_1q_all - import$MMrisk_13_all

table(import$MMrisk_CLASSIFIER_ALL)

import$MMrisk_AB_ALL <- ifelse(import$MMrisk_CLASSIFIER_ALL == 2 & import$MMrisk_1q_all ==1,
                               "2_1q", 
                               ifelse(import$MMrisk_CLASSIFIER_ALL == 2 & import$MMrisk_13_all ==1,
                                      "2_13",
                                      import$MMrisk_CLASSIFIER_ALL)
)





TabTxtable <- table(import$MMrisk_CLASSIFIER_ALL, import$TX_I_LINE_1_0)
TabTxtable_df <-  as.matrix.data.frame(TabTxtable)
TabTxtable_df


table(import$MMrisk_AB_ALL)

import$risk1_VS_1q_VS_other <- recode(import$MMrisk_AB_ALL, "1"="1q&13", "2_1q"="1q_only", "2_13"="other", "3"="other")
table(import$risk1_VS_1q_VS_other)

OS <- Surv(import$OS_MONTHS, import$OS_EVENT)
PFS <- Surv(import$PFS_I_MONTHS, import$PFS_I_EVENT)


#______ ULTRA MM RISK CREATION _____

import$ULTRA_MMRISK <- ifelse(import$MMrisk_CLASSIFIER_ALL==1 & (import$ISS==3 | import$ISS ==2), "Ultra_High_Risk", 
                              ifelse(import$MMrisk_CLASSIFIER_ALL==3 & import$ISS==1, "Ultra_Low_Risk", "other"))


import$ULTRA_MMRISK %>% table


#==================== SURVIVAL =====================

#_______ PROTOCOL ________ **

ggsurvplot(survfit(OS ~ import$PROTOCOL_REV, data = import), pval = T, risk.table = T)+ xlab("OS")
ggsurvplot(survfit(PFS ~ import$PROTOCOL_REV, data = import), pval = T, risk.table = T) + xlab("PFS")

#_______ OLD ________ ***

ggsurvplot(survfit(OS ~ import$OLD_0_1, data = import), pval = T, risk.table = T)+ xlab("OS")
ggsurvplot(survfit(PFS ~ import$OLD_0_1, data = import), pval = T, risk.table = T)+ xlab("PFS")

#_______ SEX ________ /

ggsurvplot(survfit(OS ~ import$SEX, data = import), pval = T, risk.table = T)+ xlab("OS")
ggsurvplot(survfit(PFS ~ import$SEX, data = import), pval = T, risk.table = T)+ xlab("PFS")


#_______ HB _______ ***

ggsurvplot(survfit(OS ~ import$HB_m_105, data = import), pval = T, risk.table = T)+ xlab("OS")
ggsurvplot(survfit(PFS ~ import$HB_m_105, data = import), pval = T, risk.table = T)+ xlab("PFS")


#_______ PLT ________ ***

ggsurvplot(survfit(OS ~ import$PLT_m_150, data = import), pval = T, risk.table = T) + xlab("OS")
ggsurvplot(survfit(PFS ~ import$PLT_m_150, data = import), pval = T, risk.table = T)+ xlab("PFS")

#________ PC ________ /

ggsurvplot(survfit(OS ~ import$PC_M_60, data = import), pval = T, risk.table = T)+ xlab("OS")
ggsurvplot(survfit(PFS ~ import$PC_M_60, data = import), pval = T, risk.table = T)+ xlab("PFS")

#______ LDH __________ /

ggsurvplot(survfit(OS ~ import$LDH_UL, data = import), pval = T, risk.table = T)+ xlab("OS")
ggsurvplot(survfit(PFS ~ import$LDH_UL, data = import), pval = T, risk.table = T)+ xlab("PFS")

#______ ISS __________ ***

ggsurvplot(survfit(OS ~ import$ISS, data = import), pval = T, risk.table = T)+ xlab("OS")
ggsurvplot(survfit(PFS ~ import$ISS, data = import), pval = T, risk.table = T)+ xlab("PFS")

#______ R-ISS __________ ***

ggsurvplot(survfit(OS ~ import$R_ISS, data = import), pval = T, risk.table = T)+ xlab("OS")
ggsurvplot(survfit(PFS ~ import$R_ISS, data = import), pval = T, risk.table = T)+ xlab("PFS")

#______ Ig_Isotype __________ /

ggsurvplot(survfit(OS ~ import$IG_ISOTYPE_REV, data = import), pval = T, risk.table = T)+ xlab("OS")
ggsurvplot(survfit(PFS ~ import$IG_ISOTYPE_REV, data = import), pval = T, risk.table = T)+ xlab("PFS")

#______ Light-Chain __________ /

ggsurvplot(survfit(OS ~ import$LIGHT_CHAIN, data = import), pval = T, risk.table = T)+ xlab("OS")
ggsurvplot(survfit(PFS ~ import$LIGHT_CHAIN, data = import), pval = T, risk.table = T)+ xlab("PFS")


#______ TX yes/no ________ ***

ggsurvplot(survfit(OS ~ import$TX_I_LINE_1_0, data = import), pval = T, risk.table = T)+ xlab("OS")
ggsurvplot(survfit(PFS ~ import$TX_I_LINE_1_0, data = import), pval = T, risk.table = T)+ xlab("PFS")

#______ number TX __________ ***

ggsurvplot(survfit(OS ~ import$NUMBERS_TX_FRONT_LINE, data = import), pval = T, risk.table = T)+ xlab("OS")
ggsurvplot(survfit(PFS ~ import$NUMBERS_TX_FRONT_LINE, data = import), pval = T, risk.table = T)+ xlab("PFS")


#_______ t(4;14) _______ ***

ggsurvplot(survfit(OS ~ import$FISH_T_4_14, data = import), pval = T, risk.table = T)+ xlab("OS")
ggsurvplot(survfit(PFS ~ import$FISH_T_4_14, data = import), pval = T, risk.table = T)+ xlab("PFS")

#_______ t(11;14) _______ /

ggsurvplot(survfit(OS ~ import$FISH_T_11_14, data = import), pval = T, risk.table = T)+ xlab("OS")
ggsurvplot(survfit(PFS ~ import$FISH_T_11_14, data = import), pval = T, risk.table = T)+ xlab("PFS")

#_______ t(14;16) _______ *

ggsurvplot(survfit(OS ~ import$FISH_T_14_16, data = import), pval = T, risk.table = T)+ xlab("OS")
ggsurvplot(survfit(PFS ~ import$FISH_T_14_16, data = import), pval = T, risk.table = T)+ xlab("PFS")

#_______ t(14;20) _______ ***

ggsurvplot(survfit(OS ~ import$FISH_T_14_20, data = import), pval = T, risk.table = T)+ xlab("OS")
ggsurvplot(survfit(PFS ~ import$FISH_T_14_20, data = import), pval = T, risk.table = T)+ xlab("PFS")

#_______ t(6;14) _______ /

ggsurvplot(survfit(OS ~ import$FISH_T_6_14, data = import), pval = T, risk.table = T)+ xlab("OS")
ggsurvplot(survfit(PFS ~ import$FISH_T_6_14, data = import), pval = T, risk.table = T)+ xlab("PFS")


#_______ Hyperdiploidy _______ **

ggsurvplot(survfit(OS ~ import$HyperDiploidy, data = import), pval = T, risk.table = T)+ xlab("OS")
ggsurvplot(survfit(PFS ~ import$HyperDiploidy, data = import), pval = T, risk.table = T)+ xlab("PFS")

#______ DEL 17p _______ ***

ggsurvplot(survfit(OS ~ import$`DEL_maj-broad_chr_17p`, data = import), pval = T, risk.table = T)+ xlab("OS")
ggsurvplot(survfit(PFS ~ import$`DEL_maj-broad_chr_17p`, data = import), pval = T, risk.table = T)+ xlab("PFS")

coxph(OS ~ import$DEL_maj-broad_chr_17p,data = import) %>% summary
coxph(PFS ~ import$DEL_maj-broad_chr_17p,data = import) %>% summary

#______ DEL TP53 ________ ***

ggsurvplot(survfit(OS ~ import$`DEL_maj-focal_TP53`, data = import), pval = T, risk.table = T)+ xlab("OS")
ggsurvplot(survfit(PFS ~ import$`DEL_maj-focal_TP53`, data = import), pval = T, risk.table = T)+ xlab("PFS")

coxph(OS ~ import$DEL_maj-focal_TP53,data = import) %>% summary
coxph(PFS ~ import$DEL_maj-focal_TP53,data = import) %>% summary

#______ AMP 1q ________ ***

# Broad
ggsurvplot(survfit(OS ~ import$`AMP_maj-broad_chr_1q`, data = import), pval = T, risk.table = T)+ xlab("OS")
ggsurvplot(survfit(PFS ~ import$`AMP_maj-broad_chr_1q`, data = import), pval = T, risk.table = T)+ xlab("PFS")

coxph(OS ~ import$AMP_maj-broad_chr_1q,data = import) %>% summary
coxph(PFS ~ import$AMP_maj-broad_chr_1q,data = import) %>% summary

# 1q MMrisk focal+broad 

ggsurvplot(survfit(OS ~ import$MMrisk_1q_all, data = import), pval = T, risk.table = T)+ xlab("OS")
ggsurvplot(survfit(PFS ~ import$MMrisk_1q_all, data = import), pval = T, risk.table = T)+ xlab("PFS")

coxph(OS ~ import$MMrisk_1q_all,data = import) %>% summary
coxph(PFS ~ import$MMrisk_1q_all,data = import) %>% summary

#______ DEL 13 ________ ***

# Broad
ggsurvplot(survfit(OS ~ import$`DEL_maj-broad_chr_13q`, data = import), pval = T, risk.table = T) + xlab("OS")
ggsurvplot(survfit(PFS ~ import$`DEL_maj-broad_chr_13q`, data = import), pval = T, risk.table = T)+ xlab("PFS")

coxph(OS ~ import$DEL_maj-broad_chr_13q,data = import) %>% summary
coxph(PFS ~ import$DEL_maj-broad_chr_13q,data = import) %>% summary

# 13 MMrisk focal+broad
ggsurvplot(survfit(OS ~ import$MMrisk_13_all, data = import), pval = T, risk.table = T) + xlab("OS")
ggsurvplot(survfit(PFS ~ import$MMrisk_13_all, data = import), pval = T, risk.table = T) + xlab("PFS")

coxph(OS ~ import$MMrisk_13_all,data = import) %>% summary
coxph(PFS ~ import$MMrisk_13_all,data = import) %>% summary


#______ MM RISK ________ ***

ggsurvplot(survfit(OS ~ import$MMrisk_CLASSIFIER_ALL, data = import), pval = T, risk.table = T)+ xlab("OS")
ggsurvplot(survfit(PFS ~ import$MMrisk_CLASSIFIER_ALL, data = import), pval = T, risk.table = T)+ xlab("PFS")

coxph(OS ~ import$MMrisk_CLASSIFIER_ALL,data = import) %>% summary
coxph(PFS ~ import$MMrisk_CLASSIFIER_ALL,data = import) %>% summary


import$MMrisk1 <- ifelse(import$MMrisk_CLASSIFIER_ALL==1, 1, 0)
import$MMrisk2 <- ifelse(import$MMrisk_CLASSIFIER_ALL==2, 1, 0)
import$MMrisk3 <- ifelse(import$MMrisk_CLASSIFIER_ALL==3, 1, 0)

#______ Risk 1 ______ ***
ggsurvplot(survfit(OS ~ import$MMrisk1, data = import), pval = T, risk.table = T)+ xlab("OS")
ggsurvplot(survfit(PFS ~ import$MMrisk1, data = import), pval = T, risk.table = T)+ xlab("PFS")

coxph(OS ~ import$MMrisk1,data = import) %>% summary
coxph(PFS ~ import$MMrisk1,data = import) %>% summary

#_______ Risk 2 ______ /
ggsurvplot(survfit(OS ~ import$MMrisk2, data = import), pval = T, risk.table = T)+ xlab("OS")
ggsurvplot(survfit(PFS ~ import$MMrisk2, data = import), pval = T, risk.table = T)+ xlab("PFS")

coxph(OS ~ import$MMrisk2,data = import) %>% summary
coxph(PFS ~ import$MMrisk2,data = import) %>% summary

#_______ Risk 3 _______ ***
ggsurvplot(survfit(OS ~ import$MMrisk3, data = import), pval = T, risk.table = T)+ xlab("OS")
ggsurvplot(survfit(PFS ~ import$MMrisk3, data = import), pval = T, risk.table = T)+ xlab("PFS")

coxph(OS ~ import$MMrisk3,data = import) %>% summary
coxph(PFS ~ import$MMrisk3,data = import) %>% summary

#____ MM RISK a/b _______

ggsurvplot(survfit(OS ~ import$MMrisk_AB_ALL, data = import), pval = T, risk.table = T)+ xlab("OS")
ggsurvplot(survfit(PFS ~ import$MMrisk_AB_ALL, data = import), pval = T, risk.table = T)+ xlab("PFS")

coxph(OS ~ import$MMrisk_AB_ALL,data = import) %>% summary
coxph(PFS ~ import$MMrisk_AB_ALL,data = import) %>% summary

#____ Risk_1 vs only_1q _______

ggsurvplot(survfit(OS ~ import$risk1_VS_1q_VS_other, data = import), pval = T, risk.table = T)+ xlab("OS")
ggsurvplot(survfit(PFS ~ import$risk1_VS_1q_VS_other, data = import), pval = T, risk.table = T)+ xlab("PFS")

imp1q <- filter(import, `AMP_maj-broad_chr_1q` == 1)

OS1q <- Surv(imp1q$OS_MONTHS, imp1q$OS_EVENT)
PFS1q <- Surv(imp1q$PFS_I_MONTHS, imp1q$PFS_I_EVENT)

ggsurvplot(survfit(OS1q ~ imp1q$risk1_VS_1q_VS_other, data = imp1q), pval = T, risk.table = T)+ xlab("OS")
ggsurvplot(survfit(PFS1q ~ imp1q$risk1_VS_1q_VS_other, data = imp1q), pval = T, risk.table = T)+ xlab("PFS")

coxph(OS1q ~ imp1q$risk1_VS_1q_VS_other,data = imp1q) %>% summary
coxph(PFS1q ~ imp1q$risk1_VS_1q_VS_other,data = imp1q) %>% summary


#______ ULTRA MMrisk ________


ggsurvplot(survfit(OS ~ import$ULTRA_MMRISK, data = import), pval = T, risk.table = T)+ xlab("OS")
ggsurvplot(survfit(PFS ~ import$ULTRA_MMRISK, data = import), pval = T, risk.table = T)+ xlab("PFS")

survfit(OS ~ import$ULTRA_MMRISK)
survfit(PFS ~ import$ULTRA_MMRISK)

import$Ultra_High <- ifelse(import$ULTRA_MMRISK=="Ultra_High_Risk", 1, 0)
import$Ultra_Low <- ifelse(import$ULTRA_MMRISK=="Ultra_Low_Risk", 1, 0)


#_____ ultra High _________

ggsurvplot(survfit(OS ~ import$Ultra_High, data = import), pval = T, risk.table = T)+ xlab("OS")
ggsurvplot(survfit(PFS ~ import$Ultra_High, data = import), pval = T, risk.table = T)+ xlab("PFS")

survfit(OS ~ import$Ultra_High)
survfit(PFS ~ import$Ultra_High)

coxph(OS ~ import$Ultra_High + strata(PROTOCOL_REV), data = import) %>% summary
coxph(PFS ~ import$Ultra_High + strata(PROTOCOL_REV), data = import) %>% summary


#_____ ultra Low _________

ggsurvplot(survfit(OS ~ import$Ultra_Low, data = import), pval = T, risk.table = T)+ xlab("OS")
ggsurvplot(survfit(PFS ~ import$Ultra_Low, data = import), pval = T, risk.table = T)+ xlab("PFS")

survfit(OS ~ import$Ultra_Low)
survfit(PFS ~ import$Ultra_Low)

coxph(OS ~ import$Ultra_Low + strata(PROTOCOL_REV), data = import) %>% summary
coxph(PFS ~ import$Ultra_Low + strata(PROTOCOL_REV), data = import) %>% summary




#=======================================================================
#========================== MULTIVARIATE ===============================
#=======================================================================


#_____ multivariate ______

coxph(OS ~ MMrisk1 + `DEL_maj-focal_TP53` + FISH_T_4_14 + `DEL_maj-broad_chr_1p`, data = import) %>% summary
coxph(PFS ~ MMrisk1 + `DEL_maj-focal_TP53` + FISH_T_4_14 + `DEL_maj-broad_chr_1p`, data = import) %>% summary

coxph(OS ~ strata(PROTOCOL) + MMrisk1 + `DEL_maj-focal_TP53` + FISH_T_4_14 + `DEL_maj-broad_chr_1p`, data = import) %>% summary
coxph(PFS ~ strata(PROTOCOL) + MMrisk1 + `DEL_maj-focal_TP53` + FISH_T_4_14 + `DEL_maj-broad_chr_1p`, data = import) %>% summary

# OLD_0_1  HB_m_105  PLT_m_150  TX_I_LINE_1_0  HyperDiploidy  FISH_T_4_14  FISH_T_14_16  FISH_T_14_20


coxph(OS ~ strata(PROTOCOL) + MMrisk1 + `DEL_maj-focal_TP53` + FISH_T_4_14 + `DEL_maj-broad_chr_1p` + OLD_0_1 + HB_m_105 + PLT_m_150 + TX_I_LINE_1_0 + HyperDiploidy + FISH_T_14_16 + FISH_T_14_20 + HyperDiploidy, data = import) %>% summary
coxph(PFS ~ strata(PROTOCOL) + MMrisk1 + `DEL_maj-focal_TP53` + FISH_T_4_14 + `DEL_maj-broad_chr_1p` + OLD_0_1 + HB_m_105 + PLT_m_150 + TX_I_LINE_1_0 + HyperDiploidy + FISH_T_14_16 + FISH_T_14_20 + HyperDiploidy, data = import) %>% summary


coxph(OS ~ MMrisk3 + HyperDiploidy + t_11_14, data = import) %>% summary
coxph(PFS ~ MMrisk3 + HyperDiploidy + t_11_14, data = import) %>% summary

