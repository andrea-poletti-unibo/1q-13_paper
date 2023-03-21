library(data.table)
library(tidyverse)
library(survival)
library(survminer)

import <- fread("D:/analisi_in_corso/risk_1q_13/complete_database_1q_13_181019.txt") 
import$PROTOCOL %>% table

import$Age %>% hist(col="blue")
import$Tx_number %>% hist(col="blue")

#______ MM RISK CREATION _____

import$AMP_1q_genes_all <- with(import, 
                            ifelse(`AMP_all_focal_ANP32E`==1 | `AMP_all_focal_MCL1` ==1 | `AMP_all_focal_CKS1B`==1, 1,0 ))


table(import$`AMP_all_broad_chr_1q`, import$AMP_1q_genes_all)
table(import$`DEL_all_broad_chr_13q`, import$`DEL_all_focal_RB1`)


import$MMrisk_1q_all <- ifelse( import$`AMP_all_broad_chr_1q` == 1 | import$AMP_1q_genes_all ==1, 1, 0)
import$MMrisk_13_all <- ifelse( import$`DEL_all_broad_chr_13q` ==1 | import$`DEL_all_focal_RB1` ==1, 1, 0)

import$MMrisk_CLASSIFIER_ALL <- 3  - import$MMrisk_1q_all - import$MMrisk_13_all

table(import$MMrisk_CLASSIFIER_ALL)

import$MMrisk_AB_ALL <- ifelse(import$MMrisk_CLASSIFIER_ALL == 2 & import$MMrisk_1q_all ==1,
                          "2a", 
                          ifelse(import$MMrisk_CLASSIFIER_ALL == 2 & import$MMrisk_13_all ==1,
                                 "2b",
                                 import$MMrisk_CLASSIFIER_ALL)
                          )






table(import$MMrisk_CLASSIFIER_ALL, import$Tx_yes_no)
table(import$MMrisk_AB_ALL)

import$risk1_VS_1q_VS_other <- recode(import$MMrisk_AB_ALL, "1"="1q&13", "2a"="1q_only", "2b"="other", "3"="other")
table(import$risk1_VS_1q_VS_other)

OS <- Surv(import$OS_months, import$OS_event_death)
PFS <- Surv(import$PFS_I_months, import$PFS_I_event)


#______ TX yes/no ________

ggsurvplot(survfit(OS ~ Tx_yes_no, data = import), pval = T)
ggsurvplot(survfit(PFS ~ Tx_yes_no, data = import), pval = T)

#_______t(4;14)_______

ggsurvplot(survfit(OS ~ import$t_4_14, data = import), pval = T)
ggsurvplot(survfit(PFS ~ import$t_4_14, data = import), pval = T)

#_______t(4;14)_______

ggsurvplot(survfit(OS ~ import$HyperDiploidy, data = import), pval = T)
ggsurvplot(survfit(PFS ~ import$HyperDiploidy, data = import), pval = T)

#______ DEL 17p ________

ggsurvplot(survfit(OS ~ import$DEL_all_broad_chr_17p, data = import), pval = T)
ggsurvplot(survfit(PFS ~ import$DEL_all_broad_chr_17p, data = import), pval = T)

coxph(OS ~ import$DEL_all_broad_chr_17p,data = import) %>% summary
coxph(PFS ~ import$DEL_all_broad_chr_17p,data = import) %>% summary

#______ DEL TP53 ________

ggsurvplot(survfit(OS ~ import$DEL_all_focal_TP53, data = import), pval = T)
ggsurvplot(survfit(PFS ~ import$DEL_all_focal_TP53, data = import), pval = T)

coxph(OS ~ import$DEL_all_focal_TP53,data = import) %>% summary
coxph(PFS ~ import$DEL_all_focal_TP53,data = import) %>% summary

#______ AMP 1q ________

# Broad
ggsurvplot(survfit(OS ~ import$AMP_all_broad_chr_1q, data = import), pval = T)
ggsurvplot(survfit(PFS ~ import$AMP_all_broad_chr_1q, data = import), pval = T)

coxph(OS ~ import$AMP_all_broad_chr_1q,data = import) %>% summary
coxph(PFS ~ import$AMP_all_broad_chr_1q,data = import) %>% summary

# 1q MMrisk focal+broad 

ggsurvplot(survfit(OS ~ import$MMrisk_1q_all, data = import), pval = T)
ggsurvplot(survfit(PFS ~ import$MMrisk_1q_all, data = import), pval = T)

coxph(OS ~ import$MMrisk_1q_all,data = import) %>% summary
coxph(PFS ~ import$MMrisk_1q_all,data = import) %>% summary

#______ DEL 13 ________

# Broad
ggsurvplot(survfit(OS ~ import$DEL_all_broad_chr_13q, data = import), pval = T) # WUT
ggsurvplot(survfit(PFS ~ import$DEL_all_broad_chr_13q, data = import), pval = T) # WUT

coxph(OS ~ import$DEL_all_broad_chr_13q,data = import) %>% summary
coxph(PFS ~ import$DEL_all_broad_chr_13q,data = import) %>% summary

# 13 MMrisk focal+broad
ggsurvplot(survfit(OS ~ import$MMrisk_13_all, data = import), pval = T) # WUT
ggsurvplot(survfit(PFS ~ import$MMrisk_13_all, data = import), pval = T) # WUT

coxph(OS ~ import$MMrisk_13_all,data = import) %>% summary
coxph(PFS ~ import$MMrisk_13_all,data = import) %>% summary


#______ MM RISK ________

ggsurvplot(survfit(OS ~ import$MMrisk_CLASSIFIER_ALL, data = import), pval = T)
ggsurvplot(survfit(PFS ~ import$MMrisk_CLASSIFIER_ALL, data = import), pval = T)

coxph(OS ~ import$MMrisk_CLASSIFIER_ALL,data = import) %>% summary
coxph(PFS ~ import$MMrisk_CLASSIFIER_ALL,data = import) %>% summary


import$MMrisk1 <- ifelse(import$MMrisk_CLASSIFIER_ALL==1, 1, 0)
import$MMrisk2 <- ifelse(import$MMrisk_CLASSIFIER_ALL==2, 1, 0)
import$MMrisk3 <- ifelse(import$MMrisk_CLASSIFIER_ALL==3, 1, 0)

# Risk 1
ggsurvplot(survfit(OS ~ import$MMrisk1, data = import), pval = T)
ggsurvplot(survfit(PFS ~ import$MMrisk1, data = import), pval = T)

coxph(OS ~ import$MMrisk1,data = import) %>% summary
coxph(PFS ~ import$MMrisk1,data = import) %>% summary

# Risk 2
ggsurvplot(survfit(OS ~ import$MMrisk2, data = import), pval = T)
ggsurvplot(survfit(PFS ~ import$MMrisk2, data = import), pval = T)

coxph(OS ~ import$MMrisk2,data = import) %>% summary
coxph(PFS ~ import$MMrisk2,data = import) %>% summary

# Risk 3
ggsurvplot(survfit(OS ~ import$MMrisk3, data = import), pval = T)
ggsurvplot(survfit(PFS ~ import$MMrisk3, data = import), pval = T)

coxph(OS ~ import$MMrisk3,data = import) %>% summary
coxph(PFS ~ import$MMrisk3,data = import) %>% summary

#____ MM RISK a/b _______

ggsurvplot(survfit(OS ~ import$MMrisk_AB_ALL, data = import), pval = T)
ggsurvplot(survfit(PFS ~ import$MMrisk_AB_ALL, data = import), pval = T)

coxph(OS ~ import$MMrisk_AB_ALL,data = import) %>% summary
coxph(PFS ~ import$MMrisk_AB_ALL,data = import) %>% summary

#____ Risk_1 vs only_1q _______

ggsurvplot(survfit(OS ~ import$risk1_VS_1q_VS_other, data = import), pval = T)
ggsurvplot(survfit(PFS ~ import$risk1_VS_1q_VS_other, data = import), pval = T)

imp1q <- filter(import, `AMP_all_broad_chr_1q` == 1)

OS1q <- Surv(imp1q$OS_months, imp1q$OS_event_death)
PFS1q <- Surv(imp1q$PFS_I_months, imp1q$PFS_I_event)

ggsurvplot(survfit(OS1q ~ imp1q$risk1_VS_1q_VS_other, data = imp1q), pval = T)
ggsurvplot(survfit(PFS1q ~ imp1q$risk1_VS_1q_VS_other, data = imp1q), pval = T)

coxph(OS1q ~ imp1q$risk1_VS_1q_VS_other,data = imp1q) %>% summary


#=======================================================================

#_____ multivariate ______

coxph(OS ~ MMrisk1 + MMrisk_1q_all, data = import) %>% summary
coxph(PFS ~ MMrisk1 + MMrisk_1q_all, data = import) %>% summary


coxph(OS ~ MMrisk1 + DEL_all_focal_TP53 + t_4_14 + `DEL_all_broad_chr_1p`, data = import) %>% summary
coxph(PFS ~ MMrisk1 + DEL_all_focal_TP53 + t_4_14 + DEL_all_broad_chr_1p, data = import) %>% summary

coxph(OS ~ MMrisk3 + HyperDiploidy + t_11_14, data = import) %>% summary
coxph(PFS ~ MMrisk3 + HyperDiploidy + t_11_14, data = import) %>% summary


#______ ultra MMrisk model ______

import$ULTRA_MMrisk <- ifelse(import$ISS %in% c(3,2) & import$MMrisk_CLASSIFIER_ALL == 1,
                              1,
                              ifelse(import$ISS==1 & import$MMrisk_CLASSIFIER_ALL == 3,
                                     3,
                                     2))
import$ULTRA_MMrisk %>% table

ggsurvplot(survfit(OS ~ import$ULTRA_MMrisk, data = import), pval = T)
ggsurvplot(survfit(PFS ~ import$ULTRA_MMrisk, data = import), pval = T)



