library(data.table)
library(tidyverse)
library(survival)
library(survminer)

import <- fread("C:/Users/mm_gr/Alma Mater Studiorum Università di Bologna/PROJECT 1q & 13 - Documenti/complete_database_1q_13.txt") 
import$PROTOCOL %>% table

import$Age %>% hist(col="blue")
import$Tx_number %>% hist(col="blue")


import$MMrisk <- 3 - import$AMP_broad_chr_1q - import$DEL_broad_chr_13q

import$MMriskAB <- ifelse(import$MMrisk == 2 & import$AMP_broad_chr_1q ==1,
                          "2a", 
                          ifelse(import$MMrisk == 2 & import$DEL_broad_chr_13q ==1,
                                 "2b",
                                 import$MMrisk)
                          )

table(import$MMrisk, import$Tx_yes_no)
table(import$MMriskAB)

import$risk1_VS_1q_VS_other <- recode(import$MMriskAB, "1"="1q&13", "2a"="1q_only", "2b"="other", "3"="other")
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

ggsurvplot(survfit(OS ~ import$DEL_broad_chr_17p, data = import), pval = T)
ggsurvplot(survfit(PFS ~ import$DEL_broad_chr_17p, data = import), pval = T)

coxph(OS ~ import$DEL_broad_chr_17p,data = import) %>% summary
coxph(PFS ~ import$DEL_broad_chr_17p,data = import) %>% summary

#______ AMP 1q ________

ggsurvplot(survfit(OS ~ import$AMP_broad_chr_1q, data = import), pval = T)
ggsurvplot(survfit(PFS ~ import$AMP_broad_chr_1q, data = import), pval = T)

coxph(OS ~ import$AMP_broad_chr_1q,data = import) %>% summary
coxph(PFS ~ import$AMP_broad_chr_1q,data = import) %>% summary

#______ DEL 13 ________

ggsurvplot(survfit(OS ~ import$DEL_broad_chr_13q, data = import), pval = T) # WUT
ggsurvplot(survfit(PFS ~ import$DEL_broad_chr_13q, data = import), pval = T) # WUT

coxph(OS ~ import$DEL_broad_chr_13q,data = import) %>% summary
coxph(PFS ~ import$DEL_broad_chr_13q,data = import) %>% summary

#______ MM RISK ________

ggsurvplot(survfit(OS ~ import$MMrisk, data = import), pval = T)
ggsurvplot(survfit(PFS ~ import$MMrisk, data = import), pval = T)

coxph(OS ~ import$MMrisk,data = import) %>% summary
coxph(PFS ~ import$MMrisk,data = import) %>% summary


import$MMrisk1 <- ifelse(import$MMrisk==1, 1, 0)
import$MMrisk2 <- ifelse(import$MMrisk==2, 1, 0)
import$MMrisk3 <- ifelse(import$MMrisk==3, 1, 0)

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

ggsurvplot(survfit(OS ~ import$MMriskAB, data = import), pval = T)
ggsurvplot(survfit(PFS ~ import$MMriskAB, data = import), pval = T)

coxph(OS ~ import$MMriskAB,data = import) %>% summary
coxph(PFS ~ import$MMriskAB,data = import) %>% summary

#____ Risk_1 vs only_1q _______

ggsurvplot(survfit(OS ~ import$risk1_VS_1q_VS_other, data = import), pval = T)
ggsurvplot(survfit(PFS ~ import$risk1_VS_1q_VS_other, data = import), pval = T)

imp1q <- filter(import, AMP_broad_chr_1q == 1)

OS1q <- Surv(imp1q$OS_months, imp1q$OS_event_death)
PFS1q <- Surv(imp1q$PFS_I_months, imp1q$PFS_I_event)

ggsurvplot(survfit(OS1q ~ imp1q$risk1_VS_1q_VS_other, data = imp1q), pval = T)
ggsurvplot(survfit(PFS1q ~ imp1q$risk1_VS_1q_VS_other, data = imp1q), pval = T)

coxph(OS1q ~ imp1q$risk1_VS_1q_VS_other,data = imp1q) %>% summary


#=======================================================================

#_____ multivariate ______

coxph(OS ~ MMrisk1 + AMP_broad_chr_1q, data = import) %>% summary
coxph(PFS ~ MMrisk1 + AMP_broad_chr_1q, data = import) %>% summary


coxph(OS ~ MMrisk1 + DEL_broad_chr_17p + t_4_14 + DEL_broad_chr_1p, data = import) %>% summary
coxph(PFS ~ MMrisk1 + DEL_broad_chr_17p + t_4_14 + DEL_broad_chr_1p, data = import) %>% summary

coxph(OS ~ MMrisk3 + HyperDiploidy + t_11_14, data = import) %>% summary
coxph(PFS ~ MMrisk3 + HyperDiploidy + t_11_14, data = import) %>% summary


#______ ultra MMrisk model ______

import$ULTRA_MMrisk <- ifelse(import$ISS %in% c(3,2) & import$MMrisk == 1,
                              1,
                              ifelse(import$ISS==1 & import$MMrisk == 3,
                                     3,
                                     2))
import$ULTRA_MMrisk %>% table

ggsurvplot(survfit(OS ~ import$ULTRA_MMrisk, data = import), pval = T)
ggsurvplot(survfit(PFS ~ import$ULTRA_MMrisk, data = import), pval = T)

coxph(OS ~ MMrisk1 + AMP_broad_chr_1q, data = import) %>% summary
coxph(PFS ~ MMrisk1 + AMP_broad_chr_1q, data = import) %>% summary


