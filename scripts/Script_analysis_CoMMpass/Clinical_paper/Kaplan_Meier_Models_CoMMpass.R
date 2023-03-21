
library(tidyverse)
library(survival)
library(survminer)


#______import clinical data already created
import <- data.table::fread("clinical_data_COMMPASS-IA13_210122.txt") %>% as.data.frame()

outpath <- "C:/Users/andre/Alma Mater Studiorum Università di Bologna/PROJECT 1q & 13 - Documenti/PAPER_FIGURES_OFFICIAL/kaplan_meier/"
dir.create(outpath)

#_____ADDITIONAL DATA MANAGEMENT____

#============ create trasnlocations vars ============
import$IgH_translocation_type <-ifelse(import$SeqWGS_WHSC1_CALL==1, "t(4;14)",
                                       ifelse(import$SeqWGS_CCND1_CALL==1, "t(11;14)",
                                              ifelse(import$SeqWGS_CCND3_CALL==1,"t(6;14)",
                                                     ifelse(import$SeqWGS_MAF_CALL==1,"t(14;16)",
                                                            ifelse(import$SeqWGS_MAFB_CALL==1,"t(14;20)", 
                                                                   "no translocation")))))

import$IgH_translocation_type<-as.factor(import$IgH_translocation_type)

import$SeqFISH_T_4_14 <- ifelse(import$SeqWGS_WHSC1_CALL==1, 1, 0)
import$SeqFISH_T_11_14 <- ifelse(import$SeqWGS_CCND1_CALL==1, 1, 0)
import$SeqFISH_T_14_16 <- ifelse(import$SeqWGS_MAF_CALL==1, 1, 0)
import$SeqFISH_T_14_20 <- ifelse(import$SeqWGS_MAFB_CALL==1, 1, 0)
import$SeqFISH_T_6_14 <- ifelse(import$SeqWGS_CCND3_CALL==1, 1, 0)


#========== MM RISK CREATION ============
import$AMP_1q_genes_all <- with(import, 
                                ifelse(`AMP_maj_focal_ANP32E`==1 | `AMP_maj_focal_MCL1` ==1 | `AMP_maj_focal_CKS1B`==1, 1,0 ))

import$MMrisk_1q_all<- ifelse( import$`AMP_maj_broad_chr_1q` == 1 | import$AMP_1q_genes_all ==1, 1, 0)
import$MMrisk_13_all <- ifelse( import$`DEL_maj_broad_chr_13q` ==1 | import$`DEL_maj_focal_RB1` ==1, 1, 0)

import$MMrisk_CLASS <- 3  - import$MMrisk_1q_all - import$MMrisk_13_all
import$MMrisk_class <- paste0("risk_",import$MMrisk_CLASS) %>% as.factor() %>% relevel("risk_3")

import$MMrisk_class <- import$MMrisk_class %>% recode_factor(risk_1 ="1q&13+", risk_2 ="1q/13", risk_3 ="1q&13-") %>% relevel("1q&13-")
import$MMrisk_class %>% table()

import$MMrisk_allclass <- ifelse(import$MMrisk_class=="1q/13" & import$MMrisk_1q_all==1, "gain 1q only ", 
                                 ifelse(import$MMrisk_class=="1q/13" & import$MMrisk_13_all==1, "del 13q only",import$MMrisk_class %>% as.character())) %>% as.factor()

import$MMrisk_allclass

import$MMrisk1 <- ifelse(import$MMrisk_CLASS==1, 1, 0)
import$MMrisk2 <- ifelse(import$MMrisk_CLASS==2, 1, 0)
import$MMrisk3 <- ifelse(import$MMrisk_CLASS==3, 1, 0)

import$MMrisk_AB_ALL <- ifelse(import$MMrisk_CLASS == 2 & import$MMrisk_1q_all ==1,
                               "2_1q", 
                               ifelse(import$MMrisk_CLASS == 2 & import$MMrisk_13_all ==1,
                                      "2_13",
                                      import$MMrisk_CLASS))

import$MMrisk_AB_ALL<-as.factor(import$MMrisk_AB_ALL)
import$MMrisk_AB_ALL<-relevel(import$MMrisk_AB_ALL, ref=4)
levels(import$MMrisk_AB_ALL)


import$risk1_VS_1q_VS_other <- recode(import$MMrisk_AB_ALL, "1"="1q&13", "2_1q"="1q_only", "2_13"="other", "3"="other")


#===== create pure vars MMrisk and pure CCND2 =======

import$CCND2_traslocation <- 0
import$CCND2_traslocation[import$SeqFISH_T_4_14==1] <- 1
import$CCND2_traslocation[import$SeqFISH_T_14_16==1] <- 1
import$CCND2_traslocation[import$SeqFISH_T_14_20==1] <- 1

import$MMrisk_class_t_CCND2 <- ifelse(import$MMrisk_CLASS==1 & import$CCND2_traslocation ==1 , "t&1q&13+",
                                      ifelse(import$MMrisk_CLASS==1, "1q&13+_pure", import$MMrisk_allclass %>% as.character))

# import %>% select(MMrisk_class_t_CCND2, SeqFISH_T_4_14, SeqFISH_T_14_16, SeqFISH_T_14_20, MMrisk_CLASS) %>% View                                         

import$SeqFISH_T_4_14_pure <- ifelse(import$MMrisk_CLASS != 1 &  import$SeqFISH_T_4_14 ==1, 1, 0)
import$SeqFISH_T_14_16_pure <- ifelse(import$MMrisk_CLASS != 1 & import$SeqFISH_T_14_16 ==1, 1, 0)
import$SeqFISH_T_14_20_pure <- ifelse(import$MMrisk_CLASS != 1 & import$SeqFISH_T_14_20 ==1, 1, 0)


#========== Create Survival Data ============

import$OS_months <- import$ttcos / 30.5
import$PFS_months <- import$ttcpfs / 30.5

OS <- Surv(import$OS_months, import$censos)
PFS <- Surv(import$PFS_months, import$censpfs)


########################### KAPLAN MEIER #############################

#====================== MODEL 1: MMrisk_allclass ==== ==========================
gg <- ggsurvplot(survfit(OS ~ import$MMrisk_allclass, data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)

ggsave(plot = print(gg), filename = "KM_surv_OS_MMrisk_allclass.png", path = outpath, dpi = 300, height = 8, width = 10)

gg <- ggsurvplot(survfit(PFS ~ import$MMrisk_allclass, data = import) , pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)

ggsave(plot = print(gg), filename = "KM_surv_PFS_MMrisk_allclass.png", path = outpath, dpi = 300, height = 8, width = 10)


#====================== MODEL 2: MMrisk_class_t_CCND2 ==========================
gg <- ggsurvplot(survfit(OS ~ import$MMrisk_class_t_CCND2, data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)

ggsave(plot = print(gg), filename = "KM_surv_OS_MMrisk_class_t_CCND2.png", path = outpath, dpi = 300, height = 9, width = 12)

gg <- ggsurvplot(survfit(PFS ~ import$MMrisk_class_t_CCND2, data = import) , pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)

ggsave(plot = print(gg), filename = "KM_surv_PFS_MMrisk_class_t_CCND2.png", path = outpath, dpi = 300, height = 9, width = 12)
