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
                                ifelse(`AMP_maj_focal_ANP32E`==1 | `AMP_maj_focal_MCL1` ==1 | `AMP_maj_focal_CKS1B`==1, 1,0 ))


table(import$`AMP_maj_broad_chr_1q`, import$AMP_1q_genes_all)
table(import$`DEL_maj_broad_chr_13q`, import$`DEL_maj_focal_RB1`)


import$MMrisk_1q_all <- ifelse( import$`AMP_maj_broad_chr_1q` == 1 | import$AMP_1q_genes_all ==1, 1, 0)
import$MMrisk_13_all <- ifelse( import$`DEL_maj_broad_chr_13q` ==1 | import$`DEL_maj_focal_RB1` ==1, 1, 0)

import$MMrisk_CLASSIFIER_ALL <- 3  - import$MMrisk_1q_all - import$MMrisk_13_all

table(import$MMrisk_CLASSIFIER_ALL)

import$MMrisk_AB_ALL <- ifelse(import$MMrisk_CLASSIFIER_ALL == 2 & import$MMrisk_1q_all ==1,
                               "2a", 
                               ifelse(import$MMrisk_CLASSIFIER_ALL == 2 & import$MMrisk_13_all ==1,
                                      "2b",
                                      import$MMrisk_CLASSIFIER_ALL)
)





TabTxtable <- table(import$MMrisk_CLASSIFIER_ALL, import$Tx_yes_no)
TabTxtable_df <-  as.matrix.data.frame(TabTxtable)



table(import$MMrisk_AB_ALL)

import$risk1_VS_1q_VS_other <- recode(import$MMrisk_AB_ALL, "1"="1q&13", "2a"="1q_only", "2b"="other", "3"="other")
table(import$risk1_VS_1q_VS_other)

OS <- Surv(import$OS_months, import$OS_event_death)
PFS <- Surv(import$PFS_I_months, import$PFS_I_event)


#______ SEX _______

ggsurvplot(survfit(OS ~ Sex, data = import), pval = T, risk.table = T)
ggsurvplot(survfit(PFS ~ Sex, data = import), pval = T, risk.table = T)


#____ sex in risk groups ____
import$R3_M_F <- ifelse(import$MMrisk_CLASSIFIER_ALL==3 & import$Sex == "F", "F",
                        ifelse(import$MMrisk_CLASSIFIER_ALL==3 & import$Sex == "M", "M", NA))

import$R3_M_F %>% table

ggsurvplot(survfit(OS ~ R3_M_F, data = import), pval = T, risk.table = T)
ggsurvplot(survfit(PFS ~ R3_M_F, data = import), pval = T, risk.table = T)


# FEMALES ONLY
import$FEMALES_MMRISK <- ifelse(import$Sex=="F" & import$MMrisk_CLASSIFIER_ALL==1, "MM1_F", 
                                ifelse(import$Sex=="F" & import$MMrisk_CLASSIFIER_ALL==2, "MM2_F",
                                       ifelse(import$Sex=="F" & import$MMrisk_CLASSIFIER_ALL==3, "MM3_F", NA)))

import$FEMALES_MMRISK %>% table

ggsurvplot(survfit(OS ~ FEMALES_MMRISK, data = import), pval = T, risk.table = T)
ggsurvplot(survfit(PFS ~ FEMALES_MMRISK, data = import), pval = T, risk.table = T)


# MALES ONLY
import$MALES_MMRISK <- ifelse(import$Sex=="M" & import$MMrisk_CLASSIFIER_ALL==1, "MM1_M", 
                                ifelse(import$Sex=="M" & import$MMrisk_CLASSIFIER_ALL==2, "MM2_M",
                                       ifelse(import$Sex=="M" & import$MMrisk_CLASSIFIER_ALL==3, "MM3_M", NA)))

import$MALES_MMRISK %>% table

ggsurvplot(survfit(OS ~ MALES_MMRISK, data = import), pval = T, risk.table = T)
ggsurvplot(survfit(PFS ~ MALES_MMRISK, data = import), pval = T, risk.table = T)


#______ TX yes/no ________

ggsurvplot(survfit(OS ~ Tx_yes_no, data = import), pval = T, risk.table = T)
ggsurvplot(survfit(PFS ~ Tx_yes_no, data = import), pval = T, risk.table = T)

#_______t(4;14)_______

ggsurvplot(survfit(OS ~ import$t_4_14, data = import), pval = T, risk.table = T)
ggsurvplot(survfit(PFS ~ import$t_4_14, data = import), pval = T, risk.table = T)

#_______t(4;14)_______

ggsurvplot(survfit(OS ~ import$HyperDiploidy, data = import), pval = T, risk.table = T)
ggsurvplot(survfit(PFS ~ import$HyperDiploidy, data = import), pval = T, risk.table = T)

#______ DEL 17p ________

ggsurvplot(survfit(OS ~ import$DEL_maj_broad_chr_17p, data = import), pval = T, risk.table = T)
ggsurvplot(survfit(PFS ~ import$DEL_maj_broad_chr_17p, data = import), pval = T, risk.table = T)

coxph(OS ~ import$DEL_maj_broad_chr_17p,data = import) %>% summary
coxph(PFS ~ import$DEL_maj_broad_chr_17p,data = import) %>% summary

#______ DEL TP53 ________

ggsurvplot(survfit(OS ~ import$DEL_maj_focal_TP53, data = import), pval = T, risk.table = T)
ggsurvplot(survfit(PFS ~ import$DEL_maj_focal_TP53, data = import), pval = T, risk.table = T)

coxph(OS ~ import$DEL_maj_focal_TP53,data = import) %>% summary
coxph(PFS ~ import$DEL_maj_focal_TP53,data = import) %>% summary

#______ AMP 1q ________

# Broad
ggsurvplot(survfit(OS ~ import$AMP_maj_broad_chr_1q, data = import), pval = T, risk.table = T)
ggsurvplot(survfit(PFS ~ import$AMP_maj_broad_chr_1q, data = import), pval = T, risk.table = T)

coxph(OS ~ import$AMP_maj_broad_chr_1q,data = import) %>% summary
coxph(PFS ~ import$AMP_maj_broad_chr_1q,data = import) %>% summary

# 1q MMrisk focal+broad 

ggsurvplot(survfit(OS ~ import$MMrisk_1q_all, data = import), pval = T, risk.table = T)
ggsurvplot(survfit(PFS ~ import$MMrisk_1q_all, data = import), pval = T, risk.table = T)

coxph(OS ~ import$MMrisk_1q_all,data = import) %>% summary
coxph(PFS ~ import$MMrisk_1q_all,data = import) %>% summary

#______ DEL 13 ________

# Broad
ggsurvplot(survfit(OS ~ import$DEL_maj_broad_chr_13q, data = import), pval = T, risk.table = T) # WUT
ggsurvplot(survfit(PFS ~ import$DEL_maj_broad_chr_13q, data = import), pval = T, risk.table = T) # WUT

coxph(OS ~ import$DEL_maj_broad_chr_13q,data = import) %>% summary
coxph(PFS ~ import$DEL_maj_broad_chr_13q,data = import) %>% summary

# 13 MMrisk focal+broad
ggsurvplot(survfit(OS ~ import$MMrisk_13_all, data = import), pval = T, risk.table = T) # WUT
ggsurvplot(survfit(PFS ~ import$MMrisk_13_all, data = import), pval = T, risk.table = T) # WUT

coxph(OS ~ import$MMrisk_13_all,data = import) %>% summary
coxph(PFS ~ import$MMrisk_13_all,data = import) %>% summary


#______ MM RISK ________

ggsurvplot(survfit(OS ~ import$MMrisk_CLASSIFIER_ALL, data = import), pval = T, risk.table = T)
ggsurvplot(survfit(PFS ~ import$MMrisk_CLASSIFIER_ALL, data = import), pval = T, risk.table = T)

coxph(OS ~ import$MMrisk_CLASSIFIER_ALL,data = import) %>% summary
coxph(PFS ~ import$MMrisk_CLASSIFIER_ALL,data = import) %>% summary


import$MMrisk1 <- ifelse(import$MMrisk_CLASSIFIER_ALL==1, 1, 0)
import$MMrisk2 <- ifelse(import$MMrisk_CLASSIFIER_ALL==2, 1, 0)
import$MMrisk3 <- ifelse(import$MMrisk_CLASSIFIER_ALL==3, 1, 0)

# Risk 1
ggsurvplot(survfit(OS ~ import$MMrisk1, data = import), pval = T, risk.table = T)
ggsurvplot(survfit(PFS ~ import$MMrisk1, data = import), pval = T, risk.table = T)

coxph(OS ~ import$MMrisk1,data = import) %>% summary
coxph(PFS ~ import$MMrisk1,data = import) %>% summary

# Risk 2
ggsurvplot(survfit(OS ~ import$MMrisk2, data = import), pval = T, risk.table = T)
ggsurvplot(survfit(PFS ~ import$MMrisk2, data = import), pval = T, risk.table = T)

coxph(OS ~ import$MMrisk2,data = import) %>% summary
coxph(PFS ~ import$MMrisk2,data = import) %>% summary

# Risk 3
ggsurvplot(survfit(OS ~ import$MMrisk3, data = import), pval = T, risk.table = T)
ggsurvplot(survfit(PFS ~ import$MMrisk3, data = import), pval = T, risk.table = T)

coxph(OS ~ import$MMrisk3,data = import) %>% summary
coxph(PFS ~ import$MMrisk3,data = import) %>% summary

#____ MM RISK a/b _______

ggsurvplot(survfit(OS ~ import$MMrisk_AB_ALL, data = import), pval = T, risk.table = T)
ggsurvplot(survfit(PFS ~ import$MMrisk_AB_ALL, data = import), pval = T, risk.table = T)

coxph(OS ~ import$MMrisk_AB_ALL,data = import) %>% summary
coxph(PFS ~ import$MMrisk_AB_ALL,data = import) %>% summary

#____ Risk_1 vs only_1q _______

ggsurvplot(survfit(OS ~ import$risk1_VS_1q_VS_other, data = import), pval = T, risk.table = T)
ggsurvplot(survfit(PFS ~ import$risk1_VS_1q_VS_other, data = import), pval = T, risk.table = T)

imp1q <- filter(import, `AMP_maj_broad_chr_1q` == 1)

OS1q <- Surv(imp1q$OS_months, imp1q$OS_event_death)
PFS1q <- Surv(imp1q$PFS_I_months, imp1q$PFS_I_event)

ggsurvplot(survfit(OS1q ~ imp1q$risk1_VS_1q_VS_other, data = imp1q), pval = T, risk.table = T)
ggsurvplot(survfit(PFS1q ~ imp1q$risk1_VS_1q_VS_other, data = imp1q), pval = T, risk.table = T)

coxph(OS1q ~ imp1q$risk1_VS_1q_VS_other,data = imp1q) %>% summary


#=======================================================================

#_____ multivariate ______

coxph(OS ~ MMrisk1 + MMrisk_1q_all, data = import) %>% summary
coxph(PFS ~ MMrisk1 + MMrisk_1q_all, data = import) %>% summary


coxph(OS ~ MMrisk1 + DEL_maj_focal_TP53 + t_4_14 + `DEL_maj_broad_chr_1p`, data = import) %>% summary
coxph(PFS ~ MMrisk1 + DEL_maj_focal_TP53 + t_4_14 + DEL_maj_broad_chr_1p, data = import) %>% summary

coxph(OS ~ MMrisk3 + HyperDiploidy + t_11_14, data = import) %>% summary
coxph(PFS ~ MMrisk3 + HyperDiploidy + t_11_14, data = import) %>% summary


#______ ultra MMrisk model ______

import$ULTRA_MMrisk <- ifelse(import$ISS %in% c(3,2) & import$MMrisk_CLASSIFIER_ALL == 1,
                              1,
                              ifelse(import$ISS==1 & import$MMrisk_CLASSIFIER_ALL == 3,
                                     3,
                                     2))
import$ULTRA_MMrisk %>% table

ggsurvplot(survfit(OS ~ import$ULTRA_MMrisk, data = import), pval = T, risk.table = T)
ggsurvplot(survfit(PFS ~ import$ULTRA_MMrisk, data = import), pval = T, risk.table = T)



