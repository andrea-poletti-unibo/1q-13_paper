library(tidyverse)
library(survival)
library(survminer)
library(RODBC)
library(broom)


import <- data.table::fread("C:/Users/andre/Alma Mater Studiorum Università di Bologna/PROJECT 1q & 13 - Documenti/complete_database_1q_13_CoMMpass_240121.txt") 



# add sex and clincal data
import_clinical <- data.table::fread("C:/Users/andre/Alma Mater Studiorum Università di Bologna/PROJECT 1q & 13 - Documenti/MMRF_CoMMpass_IA13_PER_PATIENT.csv")

import_clinical$Study_Visit_iD <- paste0(import_clinical$PUBLIC_ID, "_1_BM")

import$Study_Visit_iD
import_clinical$Study_Visit_iD

import2 <- left_join(import, import_clinical,  by = "Study_Visit_iD")

import <- import2





#______ MM RISK CREATION _____
import$AMP_1q_genes_all <- with(import, 
                                ifelse(`AMP_maj_focal_ANP32E`==1 | `AMP_maj_focal_MCL1` ==1 | `AMP_maj_focal_CKS1B`==1, 1,0 ))

import$MMrisk_1q_all <- ifelse( import$`AMP_maj_broad_chr_1q` == 1 | import$AMP_1q_genes_all ==1, 1, 0)
import$MMrisk_13_all <- ifelse( import$`DEL_maj_broad_chr_13q` ==1 | import$`DEL_maj_focal_RB1` ==1, 1, 0)

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


import$ULTRA_MMrisk <- ifelse(import$R_ISS == 3 & import$MMrisk_CLASSIFIER_ALL == 1,
                              "Ultra_High",
                              ifelse(import$R_ISS==1 & import$MMrisk_CLASSIFIER_ALL == 3,
                                     "Ultra_Low",
                                     "Other"))
import$ULTRA_MMrisk %>% table

import$ULTRA_High <- ifelse(import$ULTRA_MMrisk=="Ultra_High",1,0)
import$ULTRA_Low  <- ifelse(import$ULTRA_MMrisk=="Ultra_Low",1,0)



#______ create surv ______



import$OS_months <- import$ttcos / 30.5
import$PFS_months <- import$ttcpfs / 30.5

OS <- Surv(import$OS_months, import$censos)
PFS <- Surv(import$PFS_months, import$censpfs)

#_____ other variables ______

import$SEX <- ifelse(import$D_PT_gender==1, "M", "F")

import$female_gender <- ifelse(import$SEX=="F", 1,0)

import$OLD <- ifelse(import$D_PT_age > 64, 1,0)


B2M_med <- import$D_LAB_serum_beta2_microglobulin %>% median(na.rm = T)
import$B2M_m_median <- ifelse(import$D_LAB_serum_beta2_microglobulin >B2M_med, 1,0)

albumin_med <- import$D_LAB_chem_albumin %>% median(na.rm = T)
import$ALBUMIN_m_median <- ifelse(import$D_LAB_chem_albumin>albumin_med, 1,0)


# hemoglobin in mmol/L - threshold 6.52 mmol/L equal to 10.5 g/dL  
import$HB_m_105 <- ifelse(import$D_LAB_cbc_hemoglobin < 6.52 ,1,0)

import$PLT_m_150 <- ifelse(import$D_LAB_cbc_platelet < 150 ,1,0)

import$PC_M_60 <- ifelse(import$BB_PERCENTOFPLAS > 60, 1,0)


################################################################################ 
################################ SURV ANALYSIS #################################
################################################################################ 


outpath <- "C:/Users/andre/Alma Mater Studiorum Università di Bologna/PROJECT 1q & 13 - Documenti/Survival_analysis/CoMMpass_dataset/analysis_180321/"


write_tsv(data.frame("Univariate"),paste0(outpath,"report_cox_analysis_190321.txt"), append = T)



#_______ OLD ________ 

#gg <- ggsurvplot(survfit(OS ~ import$OLD, data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("AGE <65 YEARS", "AGE >= 65 YEARS"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
##print(gg)
##ggsave(plot = #print(gg), filename = "OLD_OS.png", path = outpath, dpi = 300, height = 6, width = 8)
#gg <- ggsurvplot(survfit(PFS ~ import$OLD, data = import) , pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("AGE <65 YEARS", "AGE >= 65 YEARS"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
#print(gg)
#ggsave(plot = #print(gg), filename = "OLD_PFS.png", path = outpath, dpi = 300, height = 6, width = 8)

write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_190321.txt"), append = T)
with(import, coxph(OS ~ OLD  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_190321.txt"), append = T, col_names = T)
with(import, coxph(PFS ~ OLD  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_190321.txt"), append = T, col_names = T)


#_______ SEX ________ 

#gg <- ggsurvplot(survfit(OS ~ import$SEX, data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("SEX=F", "SEX=M"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
#print(gg)
#ggsave(plot = #print(gg), filename = "SEX_OS.png", path = outpath, dpi = 300, height = 6, width = 8)
#gg <- ggsurvplot(survfit(PFS ~ import$SEX, data = import) , pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("SEX=F", "SEX=M"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
#print(gg)
#ggsave(plot = #print(gg), filename = "SEX_PFS.png", path = outpath, dpi = 300, height = 6, width = 8)

write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_190321.txt"), append = T)
with(import, coxph(OS ~ SEX )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_190321.txt"), append = T, col_names = T)
with(import, coxph(PFS ~ SEX )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_190321.txt"), append = T, col_names = T)


#_______ B2M _______ 

#gg <- ggsurvplot(survfit(OS ~ import$B2M_m_median, data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("B2M high", "B2M low"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
#print(gg)
#ggsave(plot = #print(gg), filename = "B2M_median_OS.png", path = outpath, dpi = 300, height = 6, width = 8)
#gg <- ggsurvplot(survfit(PFS ~ import$HB_m_105, data = import) , pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("B2M high", "B2M low"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
#print(gg)
#ggsave(plot = #print(gg), filename = "B2M_median_PFS.png", path = outpath, dpi = 300, height = 6, width = 8)


write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_190321.txt"), append = T)
with(import, coxph(OS ~ D_LAB_serum_beta2_microglobulin  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_190321.txt"), append = T, col_names = T)
with(import, coxph(PFS ~ D_LAB_serum_beta2_microglobulin  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_190321.txt"), append = T, col_names = T)



#_______ ALBUMIN _______ 

#gg <- ggsurvplot(survfit(OS ~ import$ALBUMIN_m_median, data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("ALBUMIN high", "ALBUMIN low"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
#print(gg)
#ggsave(plot = #print(gg), filename = "ALBUMIN_median_OS.png", path = outpath, dpi = 300, height = 6, width = 8)
#gg <- ggsurvplot(survfit(PFS ~ import$ALBUMIN_m_median, data = import) , pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("ALBUMIN high", "ALBUMIN low"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
#print(gg)
#ggsave(plot = #print(gg), filename = "ALBUMIN_median_PFS.png", path = outpath, dpi = 300, height = 6, width = 8)


write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_190321.txt"), append = T)
with(import, coxph(OS ~ D_LAB_chem_albumin  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_190321.txt"), append = T, col_names = T)
with(import, coxph(PFS ~ D_LAB_chem_albumin )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_190321.txt"), append = T, col_names = T)




#_______ HB _______ 

#gg <- ggsurvplot(survfit(OS ~ import$HB_m_105, data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("HB>10.5 g/dL", "HB<10.5 g/dL"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
#print(gg)
#ggsave(plot = #print(gg), filename = "HB_m_105_OS.png", path = outpath, dpi = 300, height = 6, width = 8)
#gg <- ggsurvplot(survfit(PFS ~ import$HB_m_105, data = import) , pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("HB>10.5 g/dL", "HB<10.5 g/dL"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
#print(gg)
#ggsave(plot = #print(gg), filename = "HB_m_105_PFS.png", path = outpath, dpi = 300, height = 6, width = 8)


write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_190321.txt"), append = T)
with(import, coxph(OS ~ HB_m_105  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_190321.txt"), append = T, col_names = T)
with(import, coxph(PFS ~ HB_m_105  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_190321.txt"), append = T, col_names = T)


#_______ PLT ________ 

#gg <- ggsurvplot(survfit(OS ~ import$PLT_m_150, data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("PLT>150.000/mm³", "PLT<150.000/mm³"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
#print(gg)
#ggsave(plot = #print(gg), filename = "PLT_m_150_OS.png", path = outpath, dpi = 300, height = 6, width = 8)
#gg <- ggsurvplot(survfit(PFS ~ import$PLT_m_150, data = import) , pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("PLT>150.000/mm³", "PLT<150.000/mm³"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
#print(gg)
#ggsave(plot = #print(gg), filename = "PLT_m_150_PFS.png", path = outpath, dpi = 300, height = 6, width = 8)


write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_190321.txt"), append = T)
with(import, coxph(OS ~ PLT_m_150  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_190321.txt"), append = T, col_names = T)
with(import, coxph(PFS ~ PLT_m_150  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_190321.txt"), append = T, col_names = T)

#________ PC ________ 

#gg <- ggsurvplot(survfit(OS ~ import$PC_M_60, data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("PC<60%", "PC>60%"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
#print(gg)
#ggsave(plot = #print(gg), filename = "PC_M_60_OS.png", path = outpath, dpi = 300, height = 6, width = 8)
#gg <- ggsurvplot(survfit(PFS ~ import$PC_M_60, data = import) , pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("PC<60%", "PC>60%"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
#print(gg)
#ggsave(plot = #print(gg), filename = "PC_M_60_PFS.png", path = outpath, dpi = 300, height = 6, width = 8)

write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_190321.txt"), append = T)
with(import, coxph(OS ~ PC_M_60  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_190321.txt"), append = T, col_names = T)
with(import, coxph(PFS ~ PC_M_60  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_190321.txt"), append = T, col_names = T)

#______ LDH __________ 

#gg <- ggsurvplot(survfit(OS ~ import$LDH_level, data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
#print(gg)
#ggsave(plot = #print(gg), filename = "LDH_UL_OS.png", path = outpath, dpi = 300, height = 6, width = 8)
#gg <- ggsurvplot(survfit(PFS ~ import$LDH_level, data = import) , pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
#print(gg)
#ggsave(plot = #print(gg), filename = "LDH_UL_PFS.png", path = outpath, dpi = 300, height = 6, width = 8)

write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_190321.txt"), append = T)
with(import, coxph(OS ~ D_LAB_chem_ldh  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_190321.txt"), append = T, col_names = T)
with(import, coxph(PFS ~ D_LAB_chem_ldh  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_190321.txt"), append = T, col_names = T)


#______ R ISS __________ 

#gg <- ggsurvplot(survfit(OS ~ import$R_ISS, data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("ISS 1", "ISS 2", "ISS 3"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
#print(gg)
#ggsave(plot = #print(gg), filename = "ISS_OS.png", path = outpath, dpi = 300, height = 6, width = 8)
#gg <- ggsurvplot(survfit(PFS ~ import$R_ISS, data = import) , pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("ISS 1", "ISS 2", "ISS 3"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
#print(gg)
#ggsave(plot = #print(gg), filename = "ISS_PFS.png", path = outpath, dpi = 300, height = 6, width = 8)

write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_190321.txt"), append = T)
with(import, coxph(OS ~ R_ISS  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_190321.txt"), append = T, col_names = T)
with(import, coxph(PFS ~ R_ISS  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_190321.txt"), append = T, col_names = T)

#______ Light-Chain __________ 

#gg <- ggsurvplot(survfit(OS ~ import$D_IM_LIGHT_CHAIN_BY_FLOW, data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="",  tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
#print(gg)
#ggsave(plot = #print(gg), filename = "LIGHT_CHAIN_OS.png", path = outpath, dpi = 300, height = 6, width = 8)
#gg <- ggsurvplot(survfit(PFS ~ import$D_IM_LIGHT_CHAIN_BY_FLOW, data = import) , pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="",  tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
#print(gg)
#ggsave(plot = #print(gg), filename = "LIGHT_CHAIN_PFS.png", path = outpath, dpi = 300, height = 6, width = 8)


write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_190321.txt"), append = T)
with(import, coxph(OS ~ D_IM_LIGHT_CHAIN_BY_FLOW  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_190321.txt"), append = T, col_names = T)
with(import, coxph(PFS ~ D_IM_LIGHT_CHAIN_BY_FLOW  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_190321.txt"), append = T, col_names = T)


#______ TX yes/no ________ 

#gg <- ggsurvplot(survfit(OS ~ import$sctflag, data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("TX_I_LINE= 0", "TX_I_LINE= 1"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
#print(gg)
#ggsave(plot = #print(gg), filename = "TX_I_LINE_1_0_OS.png", path = outpath, dpi = 300, height = 6, width = 8)
#gg <- ggsurvplot(survfit(PFS ~ import$sctflag, data = import) , pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("TX_I_LINE= 0", "TX_I_LINE= 1"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
#print(gg)
#ggsave(plot = #print(gg), filename = "TX_I_LINE_1_0_PFS.png", path = outpath, dpi = 300, height = 6, width = 8)

write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_190321.txt"), append = T)
with(import, coxph(OS ~ sctflag  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_190321.txt"), append = T, col_names = T)
with(import, coxph(PFS ~ sctflag  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_190321.txt"), append = T, col_names = T)

#_______ t(4;14) _______ 

#gg <- ggsurvplot(survfit(OS ~ import$SeqWGS_WHSC1_CALL, data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("FISH_T_4_14=0", "FISH_T_4_14=1" ), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
#print(gg)
#ggsave(plot = #print(gg), filename = "FISH_T_4_14_OS.png", path = outpath, dpi = 300, height = 6, width = 8)
#gg <- ggsurvplot(survfit(PFS ~ import$SeqWGS_WHSC1_CALL, data = import) , pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("FISH_T_4_14=0", "FISH_T_4_14=1"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
#print(gg)
#ggsave(plot = #print(gg), filename = "FISH_T_4_14_PFS.png", path = outpath, dpi = 300, height = 6, width = 8)

write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_190321.txt"), append = T)
with(import, coxph(OS ~ SeqWGS_WHSC1_CALL  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_190321.txt"), append = T, col_names = T)
with(import, coxph(PFS ~ SeqWGS_WHSC1_CALL  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_190321.txt"), append = T, col_names = T)

#_______ t(11;14) _______ 

#gg <- ggsurvplot(survfit(OS ~ import$SeqWGS_CCND1_CALL, data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("FISH_T_11_14=0", "FISH_T_11_14=1" ), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
#print(gg)
#ggsave(plot = #print(gg), filename = "FISH_T_11_14_OS.png", path = outpath, dpi = 300, height = 6, width = 8)
#gg <- ggsurvplot(survfit(PFS ~ import$SeqWGS_CCND1_CALL, data = import) , pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("FISH_T_11_14=0", "FISH_T_11_14=1"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
#print(gg)
#ggsave(plot = #print(gg), filename = "FISH_T_11_14_PFS.png", path = outpath, dpi = 300, height = 6, width = 8)

write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_190321.txt"), append = T)
with(import, coxph(OS ~ SeqWGS_CCND1_CALL  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_190321.txt"), append = T, col_names = T)
with(import, coxph(PFS ~ SeqWGS_CCND1_CALL  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_190321.txt"), append = T, col_names = T)

#_______ t(14;16) _______ 

#gg <- ggsurvplot(survfit(OS ~ import$SeqWGS_MAF_CALL, data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("FISH_T_14_16=0", "FISH_T_14_16=1" ), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
#print(gg)
#ggsave(plot = #print(gg), filename = "FISH_T_14_16_OS.png", path = outpath, dpi = 300, height = 6, width = 8)
#gg <- ggsurvplot(survfit(PFS ~ import$SeqWGS_MAF_CALL, data = import) , pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("FISH_T_14_16=0", "FISH_T_14_16=1"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
#print(gg)
#ggsave(plot = #print(gg), filename = "FISH_T_14_16_PFS.png", path = outpath, dpi = 300, height = 6, width = 8)

write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_190321.txt"), append = T)
with(import, coxph(OS ~ SeqWGS_MAF_CALL  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_190321.txt"), append = T, col_names = T)
with(import, coxph(PFS ~ SeqWGS_MAF_CALL  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_190321.txt"), append = T, col_names = T)

#_______ t(14;20) _______ 

#gg <- ggsurvplot(survfit(OS ~ import$SeqWGS_MAFB_CALL, data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("FISH_T_14_20=0", "FISH_T_14_20=1" ), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
#print(gg)
#ggsave(plot = #print(gg), filename = "FISH_T_14_20_OS.png", path = outpath, dpi = 300, height = 6, width = 8)
#gg <- ggsurvplot(survfit(PFS ~ import$SeqWGS_MAFB_CALL, data = import) , pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("FISH_T_14_20=0", "FISH_T_14_20=1"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
#print(gg)
#ggsave(plot = #print(gg), filename = "FISH_T_14_20_PFS.png", path = outpath, dpi = 300, height = 6, width = 8)

write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_190321.txt"), append = T)
with(import, coxph(OS ~ SeqWGS_MAFB_CALL  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_190321.txt"), append = T, col_names = T)
with(import, coxph(PFS ~ SeqWGS_MAFB_CALL  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_190321.txt"), append = T, col_names = T)

#_______ t(6;14) _______ 

#gg <- ggsurvplot(survfit(OS ~ import$SeqWGS_CCND3_CALL, data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("FISH_T_6_14=0", "FISH_T_6_14=1" ), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
#print(gg)
#ggsave(plot = #print(gg), filename = "FISH_T_6_14_OS.png", path = outpath, dpi = 300, height = 6, width = 8)
#gg <- ggsurvplot(survfit(PFS ~ import$SeqWGS_CCND3_CALL, data = import) , pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("FISH_T_6_14=0", "FISH_T_6_14=1"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
#print(gg)
#ggsave(plot = #print(gg), filename = "FISH_T_6_14_PFS.png", path = outpath, dpi = 300, height = 6, width = 8)

write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_190321.txt"), append = T)
with(import, coxph(OS ~ SeqWGS_CCND3_CALL  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_190321.txt"), append = T, col_names = T)
with(import, coxph(PFS ~ SeqWGS_CCND3_CALL  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_190321.txt"), append = T, col_names = T)

#_______ Hyperdiploidy _______ 

#gg <- ggsurvplot(survfit(OS ~ import$Hyperdiploid, data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("HyperDiploidy=0", "HyperDiploidy=1" ), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
#print(gg)
#ggsave(plot = #print(gg), filename = "HyperDiploidy_OS.png", path = outpath, dpi = 300, height = 6, width = 8)
#gg <- ggsurvplot(survfit(PFS ~ import$Hyperdiploid, data = import) , pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("HyperDiploidy=0", "HyperDiploidy=1"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
#print(gg)
#ggsave(plot = #print(gg), filename = "HyperDiploidy_PFS.png", path = outpath, dpi = 300, height = 6, width = 8)


write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_190321.txt"), append = T)
with(import, coxph(OS ~ Hyperdiploid  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_190321.txt"), append = T, col_names = T)
with(import, coxph(PFS ~ Hyperdiploid  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_190321.txt"), append = T, col_names = T)


#______ DEL 17p _______ 

#gg <- ggsurvplot(survfit(OS ~ import$"DEL_maj_broad_chr_17p", data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("Chr 17p Del= 0", "Chr 17p Del= 1" ), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
#print(gg)
#ggsave(plot = #print(gg), filename = "DEL_maj_broad_chr_17p_OS.png", path = outpath, dpi = 300, height = 6, width = 8)
#gg <- ggsurvplot(survfit(PFS ~ import$"DEL_maj_broad_chr_17p", data = import) , pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("Chr 17p Del= 0", "Chr 17p Del= 1"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
#print(gg)
#ggsave(plot = #print(gg), filename = "DEL_maj_broad_chr_17p_PFS.png", path = outpath, dpi = 300, height = 6, width = 8)

write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_190321.txt"), append = T)
with(import, coxph(OS ~ `DEL_maj_broad_chr_17p`  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_190321.txt"), append = T, col_names = T)
with(import, coxph(PFS ~ `DEL_maj_broad_chr_17p`  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_190321.txt"), append = T, col_names = T)


#______ DEL 17p _______ 

#gg <- ggsurvplot(survfit(OS ~ import$"DEL_maj_broad_chr_1p", data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("Chr 1p Del= 0", "Chr 1p Del= 1" ), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
#print(gg)
#ggsave(plot = #print(gg), filename = "DEL_maj_broad_chr_1p_OS.png", path = outpath, dpi = 300, height = 6, width = 8)
#gg <- ggsurvplot(survfit(PFS ~ import$"DEL_maj_broad_chr_1p", data = import) , pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("Chr 1p Del= 0", "Chr 1p Del= 1"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
#print(gg)
#ggsave(plot = #print(gg), filename = "DEL_maj_broad_chr_1p_PFS.png", path = outpath, dpi = 300, height = 6, width = 8)

write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_190321.txt"), append = T)
with(import, coxph(OS ~ `DEL_maj_broad_chr_1p`  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_190321.txt"), append = T, col_names = T)
with(import, coxph(PFS ~ `DEL_maj_broad_chr_1p`  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_190321.txt"), append = T, col_names = T)



#______ DEL TP53 ________ 

#gg <- ggsurvplot(survfit(OS ~ import$"DEL_maj_focal_TP53", data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("TP53 Del= 0", "TP53 Del= 1" ), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
#print(gg)
#ggsave(plot = #print(gg), filename = "DEL_maj_focal_TP53_OS.png", path = outpath, dpi = 300, height = 6, width = 8)
#gg <- ggsurvplot(survfit(PFS ~ import$"DEL_maj_focal_TP53", data = import) , pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("TP53 Del= 0", "TP53 Del= 1"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
#print(gg)
#ggsave(plot = #print(gg), filename = "DEL_maj_focal_TP53_PFS.png", path = outpath, dpi = 300, height = 6, width = 8)

write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_190321.txt"), append = T)
with(import, coxph(OS ~ `DEL_maj_focal_TP53`  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_190321.txt"), append = T, col_names = T)
with(import, coxph(PFS ~ `DEL_maj_focal_TP53`  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_190321.txt"), append = T, col_names = T)

#______ AMP 1q ________ 

#gg <- ggsurvplot(survfit(OS ~ import$"AMP_maj_broad_chr_1q", data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("Chr 1q Broad Amp= 0", "Chr 1q Broad Amp= 1" ), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
#print(gg)
#ggsave(plot = #print(gg), filename = "AMP_maj_broad_chr_1_OS.png", path = outpath, dpi = 300, height = 6, width = 8)
#gg <- ggsurvplot(survfit(PFS ~ import$"AMP_maj_broad_chr_1q", data = import) , pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("Chr 1q Broad Amp= 0", "Chr 1q Broad Amp= 1"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
#print(gg)
#ggsave(plot = #print(gg), filename = "AMP_maj_broad_chr_1_PFS.png", path = outpath, dpi = 300, height = 6, width = 8)

write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_190321.txt"), append = T)
with(import, coxph(OS ~ `AMP_maj_broad_chr_1q`  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_190321.txt"), append = T, col_names = T)
with(import, coxph(PFS ~ `AMP_maj_broad_chr_1q`  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_190321.txt"), append = T, col_names = T)



# 1q MMrisk focal+broad 

#gg <- ggsurvplot(survfit(OS ~ import$MMrisk_1q_all, data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("Chr1q Broad&Focal Amp=0", "Chr1q Broad&Focal Amp=1" ), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
#print(gg)
#ggsave(plot = #print(gg), filename = "MMrisk_1q_all_OS.png", path = outpath, dpi = 300, height = 6, width = 8)
#gg <- ggsurvplot(survfit(PFS ~ import$MMrisk_1q_all, data = import) , pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("Chr1q Broad&Focal Amp=0", "Chr1q Broad&Focal Amp=1"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
#print(gg)
#ggsave(plot = #print(gg), filename = "MMrisk_1q_all_PFS.png", path = outpath, dpi = 300, height = 6, width = 8)

write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_190321.txt"), append = T)
with(import, coxph(OS ~ MMrisk_1q_all  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_190321.txt"), append = T, col_names = T)
with(import, coxph(PFS ~ MMrisk_1q_all  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_190321.txt"), append = T, col_names = T)


#______ DEL 13 ________ 

#gg <- ggsurvplot(survfit(OS ~ import$"DEL_maj_broad_chr_13q", data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("Chr 13q Broad Del= 0", "Chr 13q Broad Del= 1" ), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
#print(gg)
#ggsave(plot = #print(gg), filename = "DEL_maj_broad_chr_13q_OS.png", path = outpath, dpi = 300, height = 6, width = 8)
#gg <- ggsurvplot(survfit(PFS ~ import$"DEL_maj_broad_chr_13q", data = import) , pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("Chr 13q Broad Del= 0", "Chr 13q Broad Del= 1"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
#print(gg)
#ggsave(plot = #print(gg), filename = "DEL_maj_broad_chr_13q_PFS.png", path = outpath, dpi = 300, height = 6, width = 8)


write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_190321.txt"), append = T)
with(import, coxph(OS ~ `DEL_maj_broad_chr_13q`  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_190321.txt"), append = T, col_names = T)
with(import, coxph(PFS ~ `DEL_maj_broad_chr_13q`  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_190321.txt"), append = T, col_names = T)


# 13 MMrisk focal+broad

#gg <- ggsurvplot(survfit(OS ~ import$MMrisk_13_all, data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("Chr 13q Broad&Focal Del= 0", "Chr 13q Broad&Focal Del= 1" ), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
#print(gg)
#ggsave(plot = #print(gg), filename = "MMrisk_13_all_OS.png", path = outpath, dpi = 300, height = 6, width = 8)
#gg <- ggsurvplot(survfit(PFS ~ import$MMrisk_13_all, data = import) , pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("Chr 13q Broad&Focal Del= 0", "Chr 13q Broad&Focal Del= 1"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
#print(gg)
#ggsave(plot = #print(gg), filename = "MMrisk_13_all_PFS.png", path = outpath, dpi = 300, height = 6, width = 8)

write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_190321.txt"), append = T)
with(import, coxph(OS ~ MMrisk_13_all  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_190321.txt"), append = T, col_names = T)
with(import, coxph(PFS ~ MMrisk_13_all  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_190321.txt"), append = T, col_names = T)

#______ MM RISK ________ 

gg <- ggsurvplot(survfit(OS ~ import$MMrisk_CLASSIFIER_ALL, data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("1q&13+", "1q/13", "1q&13-" ), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
#ggsave(plot = #print(gg), filename = "MMrisk_CLASSIFIER_ALL_OS.png", path = outpath, dpi = 300, height = 6, width = 8)
gg <- ggsurvplot(survfit(PFS ~ import$MMrisk_CLASSIFIER_ALL, data = import) , pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("1q&13+", "1q/13","1q&13-" ), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
#ggsave(plot = #print(gg), filename = "MMrisk_CLASSIFIER_ALL_PFS.png", path = outpath, dpi = 300, height = 6, width = 8)


write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_190321.txt"), append = T)
with(import, coxph(OS ~ MMrisk_CLASSIFIER_ALL  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_190321.txt"), append = T, col_names = T)
with(import, coxph(PFS ~ MMrisk_CLASSIFIER_ALL  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_190321.txt"), append = T, col_names = T)

survfit(OS ~ MMrisk_CLASSIFIER_ALL, data=import  ) %>% summary
survfit(PFS ~ MMrisk_CLASSIFIER_ALL, data=import  ) %>% summary



#______ Risk 1 ______ 

#gg <- ggsurvplot(survfit(OS ~ import$MMrisk1, data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("1q&13+ =0", "1q&13+ =1" ), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
#print(gg)
#ggsave(plot = #print(gg), filename = "MMrisk1_OS.png", path = outpath, dpi = 300, height = 6, width = 8)
#gg <- ggsurvplot(survfit(PFS ~ import$MMrisk1, data = import) , pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("1q&13+ =0", "1q&13+ =1"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
#print(gg)
#ggsave(plot = #print(gg), filename = "MMrisk1_PFS.png", path = outpath, dpi = 300, height = 6, width = 8)

write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_190321.txt"), append = T)
with(import, coxph(OS ~ MMrisk1  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_190321.txt"), append = T, col_names = T)
with(import, coxph(PFS ~ MMrisk1  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_190321.txt"), append = T, col_names = T)

#_______ Risk 2 ______ 

#gg <- ggsurvplot(survfit(OS ~ import$MMrisk2, data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("1q/13+ =0", "1q/13+ =1" ), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
#print(gg)
#ggsave(plot = #print(gg), filename = "MMrisk2_OS.png", path = outpath, dpi = 300, height = 6, width = 8)
#gg <- ggsurvplot(survfit(PFS ~ import$MMrisk2, data = import) , pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("1q/13+ =0", "1q/13+ =1"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
#print(gg)
#ggsave(plot = #print(gg), filename = "MMrisk2_PFS.png", path = outpath, dpi = 300, height = 6, width = 8)


write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_190321.txt"), append = T)
with(import, coxph(OS ~ MMrisk2  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_190321.txt"), append = T, col_names = T)
with(import, coxph(PFS ~ MMrisk2  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_190321.txt"), append = T, col_names = T)


#_______ Risk 3 _______ 

#gg <- ggsurvplot(survfit(OS ~ import$MMrisk3, data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("1q&13- =0", "1q&13- =1" ), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
#print(gg)
#ggsave(plot = #print(gg), filename = "MMrisk3_OS.png", path = outpath, dpi = 300, height = 6, width = 8)
#gg <- ggsurvplot(survfit(PFS ~ import$MMrisk3, data = import) , pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("1q&13- =0", "1q&13- =1"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
#print(gg)
#ggsave(plot = #print(gg), filename = "MMrisk3_PFS.png", path = outpath, dpi = 300, height = 6, width = 8)


write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_190321.txt"), append = T)
with(import, coxph(OS ~ MMrisk3  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_190321.txt"), append = T, col_names = T)
with(import, coxph(PFS ~ MMrisk3  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_190321.txt"), append = T, col_names = T)


#____ MM RISK a/b _______

#gg <- ggsurvplot(survfit(OS ~ import$MMrisk_AB_ALL, data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("1q&13+", "only 13 Del", "only 1q Amp", "1q&13-" ), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
#print(gg)
#ggsave(plot = #print(gg), filename = "MMrisk_AB_ALL_OS.png", path = outpath, dpi = 300, height = 6, width = 8)
#gg <- ggsurvplot(survfit(PFS ~ import$MMrisk_AB_ALL, data = import) , pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("1q&13+", "only 13 Del", "only 1q Amp", "1q&13-"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
#print(gg)
#ggsave(plot = #print(gg), filename = "MMrisk_AB_ALL_PFS.png", path = outpath, dpi = 300, height = 6, width = 8)



write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_190321.txt"), append = T)
with(import, coxph(OS ~ MMrisk_AB_ALL  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_190321.txt"), append = T, col_names = T)
with(import, coxph(PFS ~ MMrisk_AB_ALL  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_190321.txt"), append = T, col_names = T)



#____ Risk_1 vs only_1q _______

#gg <- ggsurvplot(survfit(OS ~ import$risk1_VS_1q_VS_other, data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("1q&13+", "only 1q Amp", "other" ), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
#print(gg)
#ggsave(plot = #print(gg), filename = "risk1_VS_1q_VS_other_OS.png", path = outpath, dpi = 300, height = 6, width = 8)
#gg <- ggsurvplot(survfit(PFS ~ import$risk1_VS_1q_VS_other, data = import) , pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("1q&13+", "only 1q amp", "other"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
#print(gg)
#ggsave(plot = #print(gg), filename = "risk1_VS_1q_VS_other_PFS.png", path = outpath, dpi = 300, height = 6, width = 8)



write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_190321.txt"), append = T)
with(import, coxph(OS ~ risk1_VS_1q_VS_other  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_190321.txt"), append = T, col_names = T)
with(import, coxph(PFS ~ risk1_VS_1q_VS_other  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_190321.txt"), append = T, col_names = T)



imp1q <- filter(import, import$"AMP_maj_broad_chr_1q" == 1)

OS1q <- Surv(imp1q$OS_months, imp1q$censos)
PFS1q <- Surv(imp1q$PFS_months, imp1q$censpfs)

#gg <- ggsurvplot(survfit(OS1q ~ imp1q$risk1_VS_1q_VS_other, data = imp1q) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("1q&13+", "only 1q Amp" ), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
#print(gg)
#ggsave(plot = #print(gg), filename = "risk1_VS_1q_VS_other_OS1q.png", path = outpath, dpi = 300, height = 6, width = 8)
#gg <- ggsurvplot(survfit(PFS1q ~ imp1q$risk1_VS_1q_VS_other, data = imp1q) , pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("1q&13+", "only 1q Amp"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
#print(gg)
#ggsave(plot = #print(gg), filename = "risk1_VS_1q_VS_other_PFS1q.png", path = outpath, dpi = 300, height = 6, width = 8)


write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_190321.txt"), append = T)
with(imp1q, coxph(OS1q ~ risk1_VS_1q_VS_other  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_190321.txt"), append = T, col_names = T)
with(imp1q, coxph(PFS1q ~ risk1_VS_1q_VS_other  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_190321.txt"), append = T, col_names = T)


#______ ultra MMrisk classificator ______

#gg <- ggsurvplot(survfit(OS ~ import$ULTRA_MMrisk, data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("ULTRA_MMrisk=1", "ULTRA_MMrisk=2", "ULTRA_MMrisk=3" ), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
#print(gg)
#ggsave(plot = #print(gg), filename = "ULTRA_MMrisk_OS.png", path = outpath, dpi = 300, height = 6, width = 8)
#gg <- ggsurvplot(survfit(PFS ~ import$ULTRA_MMrisk, data = import), pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("ULTRA_MMrisk=1", "ULTRA_MMrisk=2", "ULTRA_MMrisk=3"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
#print(gg)
#ggsave(plot = #print(gg), filename = "ULTRA_MMrisk_PFS.png", path = outpath, dpi = 300, height = 6, width = 8)

write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_190321.txt"), append = T)
with(import, coxph(OS ~ ULTRA_MMrisk  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_190321.txt"), append = T, col_names = T)
with(import, coxph(PFS ~ ULTRA_MMrisk  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_190321.txt"), append = T, col_names = T)



#______ ultra high  ______

#gg <- ggsurvplot(survfit(OS ~ import$ULTRA_High, data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("ULTRA_High=0", "ULTRA_High=1"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
#print(gg)
#ggsave(plot = #print(gg), filename = "ULTRA_High_OS.png", path = outpath, dpi = 300, height = 6, width = 8)
#gg <- ggsurvplot(survfit(PFS ~ import$ULTRA_High, data = import), pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("ULTRA_High=0", "ULTRA_High=1"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
#print(gg)
#ggsave(plot = #print(gg), filename = "ULTRA_High_PFS.png", path = outpath, dpi = 300, height = 6, width = 8)

write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_190321.txt"), append = T)
with(import, coxph(OS ~ ULTRA_High  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_190321.txt"), append = T, col_names = T)
with(import, coxph(PFS ~ ULTRA_High  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_190321.txt"), append = T, col_names = T)



#______ ultra low  ______

#gg <- ggsurvplot(survfit(OS ~ import$ULTRA_Low, data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("ULTRA_Low=0", "ULTRA_Low=1"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
#print(gg)
#ggsave(plot = #print(gg), filename = "ULTRA_High_OS.png", path = outpath, dpi = 300, height = 6, width = 8)
#gg <- ggsurvplot(survfit(PFS ~ import$ULTRA_Low, data = import), pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("ULTRA_Low=0", "ULTRA_Low=1"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
#print(gg)
#ggsave(plot = #print(gg), filename = "ULTRA_High_PFS.png", path = outpath, dpi = 300, height = 6, width = 8)

write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_190321.txt"), append = T)
with(import, coxph(OS ~ ULTRA_Low  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_190321.txt"), append = T, col_names = T)
with(import, coxph(PFS ~ ULTRA_Low  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_190321.txt"), append = T, col_names = T)




#=======================================================================
#========================== MULTIVARIATE ===============================
#=======================================================================


.................# ggforest2 function def #......................

ggforest2 <- function (model, data = NULL, main = "Hazard ratio", cpositions = c(0.02, 
                                                                                 0.22, 0.4), fontsize = 0.7, refLabel = "reference", noDigits = 2) 
{
  conf.high <- conf.low <- estimate <- NULL
  stopifnot(inherits(model, "coxph"))
  data <- survminer:::.get_data(model, data = data)
  terms <- attr(model$terms, "dataClasses")[-1]
  coef <- as.data.frame(tidy(model, conf.int = TRUE))
  gmodel <- glance(model)
  allTerms <- lapply(seq_along(terms), function(i) {
    var <- names(terms)[i]
    if (var %in% colnames(data)) {
      if (terms[i] %in% c("factor", "character")) {
        adf <- as.data.frame(table(data[, var]))
        cbind(var = var, adf, pos = 1:nrow(adf))
      }
      else if (terms[i] == "numeric") {
        data.frame(var = var, Var1 = "", Freq = nrow(data), 
                   pos = 1)
      }
      else {
        vars = grep(paste0("^", var, "*."), coef$term, 
                    value = TRUE)
        data.frame(var = vars, Var1 = "", Freq = nrow(data), 
                   pos = seq_along(vars))
      }
    }
    else {
      message(var, "is not found in data columns, and will be skipped.")
    }
  })
  allTermsDF <- do.call(rbind, allTerms)
  colnames(allTermsDF) <- c("var", "level", "N", "pos")
  inds <- apply(allTermsDF[, 1:2], 1, paste0, collapse = "")
  rownames(coef) <- gsub(coef$term, pattern = "`", replacement = "")
  toShow <- cbind(allTermsDF, coef[inds, ])[, c("var", "level", 
                                                "N", "p.value", "estimate", "conf.low", "conf.high", 
                                                "pos")]
  toShowExp <- toShow[, 5:7]
  toShowExp[is.na(toShowExp)] <- 0
  toShowExp <- format(exp(toShowExp), digits = noDigits)
  toShowExpClean <- data.frame(toShow, pvalue = signif(toShow[, 
                                                              4], noDigits + 1), toShowExp)
  toShowExpClean$stars <- paste0(round(toShowExpClean$p.value, 
                                       noDigits + 1), " ", ifelse(toShowExpClean$p.value < 
                                                                    0.05, "*", ""), ifelse(toShowExpClean$p.value < 0.01, 
                                                                                           "*", ""), ifelse(toShowExpClean$p.value < 0.001, "*", 
                                                                                                            ""))
  toShowExpClean$ci <- paste0("(", toShowExpClean[, "conf.low.1"], 
                              " - ", toShowExpClean[, "conf.high.1"], ")")
  toShowExpClean$estimate.1[is.na(toShowExpClean$estimate)] = refLabel
  toShowExpClean$stars[which(toShowExpClean$p.value < 0.001)] = "<0.001 ***"
  toShowExpClean$stars[is.na(toShowExpClean$estimate)] = ""
  toShowExpClean$ci[is.na(toShowExpClean$estimate)] = ""
  toShowExpClean$estimate[is.na(toShowExpClean$estimate)] = 0
  toShowExpClean$var = as.character(toShowExpClean$var)
  toShowExpClean$var[duplicated(toShowExpClean$var)] = ""
  toShowExpClean$N <- paste0("(N=", toShowExpClean$N, ")")
  toShowExpClean <- toShowExpClean[nrow(toShowExpClean):1, 
  ]
  rangeb <- range(toShowExpClean$conf.low, toShowExpClean$conf.high, 
                  na.rm = TRUE)
  breaks <- axisTicks(rangeb/2, log = TRUE, nint = 7)
  rangeplot <- rangeb
  rangeplot[1] <- rangeplot[1] - diff(rangeb)
  rangeplot[2] <- rangeplot[2] + 0.15 * diff(rangeb)
  width <- diff(rangeplot)
  y_variable <- rangeplot[1] + cpositions[1] * width
  y_nlevel <- rangeplot[1] + cpositions[2] * width
  y_cistring <- rangeplot[1] + cpositions[3] * width
  y_stars <- rangeb[2]
  x_annotate <- seq_len(nrow(toShowExpClean))
  annot_size_mm <- fontsize * as.numeric(grid::convertX(unit(theme_get()$text$size, 
                                                             "pt"), "mm"))
  p <- ggplot(toShowExpClean, aes(seq_along(var), exp(estimate))) + 
    geom_rect(aes(xmin = seq_along(var) - 0.5, xmax = seq_along(var) + 
                    0.5, ymin = exp(rangeplot[1]), ymax = exp(rangeplot[2]), 
                  fill = ordered(seq_along(var)%%2 + 1))) + scale_fill_manual(values = c("#FFFFFF33", 
                                                                                         "#00000033"), guide = "none") + geom_point(pch = 15, 
                                                                                                                                    size = 4) + geom_errorbar(aes(ymin = exp(conf.low), 
                                                                                                                                                                  ymax = exp(conf.high)), width = 0.15) + geom_hline(yintercept = 1, 
                                                                                                                                                                                                                     linetype = 3) + coord_flip(ylim = exp(rangeplot)) + 
    ggtitle(main) + scale_y_log10(name = "", labels = sprintf("%g", 
                                                              breaks), expand = c(0.02, 0.02), breaks = breaks) + 
    theme_light() + theme(panel.grid.minor.y = element_blank(), 
                          panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), 
                          legend.position = "none", panel.border = element_blank(), 
                          axis.title.y = element_blank(), axis.text.y = element_blank(), 
                          axis.ticks.y = element_blank(), plot.title = element_text(hjust = 0.5)) + 
    xlab("") + annotate(geom = "text", x = x_annotate, y = exp(y_variable), 
                        label = toShowExpClean$var, fontface = "bold", hjust = 0, 
                        size = annot_size_mm) + annotate(geom = "text", x = x_annotate, 
                                                         y = exp(y_nlevel), hjust = 0, label = toShowExpClean$level, 
                                                         vjust = -0.1, size = annot_size_mm) + annotate(geom = "text", 
                                                                                                        x = x_annotate, y = exp(y_nlevel), label = toShowExpClean$N, 
                                                                                                        fontface = "italic", hjust = 0, vjust = ifelse(toShowExpClean$level == 
                                                                                                                                                         "", 0.5, 1.1), size = annot_size_mm) + annotate(geom = "text", 
                                                                                                                                                                                                         x = x_annotate, y = exp(y_cistring), label = toShowExpClean$estimate.1, 
                                                                                                                                                                                                         size = annot_size_mm, vjust = ifelse(toShowExpClean$estimate.1 == 
                                                                                                                                                                                                                                                "reference", 0.5, -0.1)) + annotate(geom = "text", 
                                                                                                                                                                                                                                                                                    x = x_annotate, y = exp(y_cistring), label = toShowExpClean$ci, 
                                                                                                                                                                                                                                                                                    size = annot_size_mm, vjust = 1.1, fontface = "italic") + 
    annotate(geom = "text", x = x_annotate, y = exp(y_stars), 
             label = toShowExpClean$stars, size = annot_size_mm, 
             hjust = -0.2, fontface = "italic") + annotate(geom = "text", 
                                                           x = 0.5, y = exp(y_variable), label = paste0("# Events: ", 
                                                                                                        gmodel$nevent, "; Global p-value (Log-Rank): ", 
                                                                                                        format.pval(gmodel$p.value.log, eps = ".001"), " \nAIC: ", 
                                                                                                        round(gmodel$AIC, 2), "; Concordance Index: ", round(gmodel$concordance, 
                                                                                                                                                             2)), size = annot_size_mm, hjust = 0, vjust = 1.2, 
                                                           fontface = "italic")
  gt <- ggplot_gtable(ggplot_build(p))
  gt$layout$clip[gt$layout$name == "panel"] <- "off"
  ggpubr::as_ggplot(gt)
}

#..........................................................


# MV1 OS MMrisk1
mv1 <- coxph(OS ~ OLD + female_gender + D_LAB_serum_beta2_microglobulin + D_LAB_chem_albumin + HB_m_105 + PLT_m_150 +
               D_LAB_chem_ldh + MMrisk1 + sctflag + `DEL_maj_focal_TP53` + SeqWGS_MAFB_CALL , data = import )
mv1 %>% summary()

mv1 <- coxph(OS ~ female_gender + D_LAB_serum_beta2_microglobulin + PLT_m_150 +
               D_LAB_chem_ldh + MMrisk1 + sctflag, data = import)
mv1 %>% summary()


ggforest2(mv1) 
ggsave(paste0(outpath,"FOREST_plot_mmrisk1_OS.png"), main = "OS 1q&13+", width = 8, height = 8)


#MV2 PFS MMrisk1
mv2 <-  with(import, coxph(PFS ~ OLD + female_gender + D_LAB_serum_beta2_microglobulin + D_LAB_chem_albumin + HB_m_105 + PLT_m_150 +
                             D_LAB_chem_ldh + MMrisk1 + sctflag, data = import )) 
mv2 %>% summary()
mv2 <- coxph(PFS ~  female_gender + D_LAB_serum_beta2_microglobulin + HB_m_105 + D_LAB_chem_ldh + MMrisk1 + sctflag, data = import )
mv2 %>% summary()

ggforest2(mv2) 
ggsave(paste0(outpath,"FOREST_plot_mmrisk1_PFS.png"), main = "PFS 1q&13+", width = 8, height = 8)


#MV3 OS UltraHigh
mv3 <- coxph(OS ~ OLD + female_gender + HB_m_105 + PLT_m_150 +
             D_LAB_chem_ldh + ULTRA_High + sctflag + `DEL_maj_focal_TP53` + SeqWGS_MAFB_CALL , data = import  )
mv3 %>% summary()
mv3 <- coxph(OS ~ female_gender + PLT_m_150 + D_LAB_chem_ldh + ULTRA_High + sctflag, data = import)
mv3 %>% summary()

ggforest2(mv3) 
ggsave(paste0(outpath,"FOREST_plot_UltraHigh_OS.png"), main = "OS Ultra High Risk", width = 8, height = 8)


#MV4 PFS UltraHigh
mv4 <- with(import, coxph(PFS ~ OLD + female_gender + HB_m_105 + PLT_m_150 +
                            D_LAB_chem_ldh + ULTRA_High + sctflag, data = import  )) 
mv4 %>% summary()
mv4 <- coxph(PFS ~ OLD + female_gender + HB_m_105 + PLT_m_150 +
               D_LAB_chem_ldh + ULTRA_High + sctflag, data = import )
mv4 %>% summary()

ggforest2(mv4) 
ggsave(paste0(outpath,"FOREST_plot_UltraHigh_PFS.png"), main = "PFS Ultra High Risk", width = 8, height = 8)




#MV5 OS MMrisk3
mv5 <- coxph(OS ~ OLD + female_gender + D_LAB_serum_beta2_microglobulin + D_LAB_chem_albumin + HB_m_105 + PLT_m_150 +
               D_LAB_chem_ldh + MMrisk3 + sctflag + `DEL_maj_focal_TP53` + SeqWGS_MAFB_CALL , data = import  )
mv5 %>% summary()

mv5 <- coxph(OS ~ female_gender + D_LAB_serum_beta2_microglobulin + PLT_m_150 +
               D_LAB_chem_ldh + MMrisk3 + sctflag , data = import )
mv5 %>% summary()


ggforest2(mv5) 
ggsave(paste0(outpath,"FOREST_plot_mmrisk3_OS.png"), main = "OS 1q&13-", width = 8, height = 8)


#MV6 PFS MMrisk3
mv6 <- coxph(PFS ~ OLD + female_gender + D_LAB_serum_beta2_microglobulin + D_LAB_chem_albumin + HB_m_105 + PLT_m_150 +
               D_LAB_chem_ldh + MMrisk3 + sctflag, data = import )
mv6 %>% summary()

mv6 <- coxph(PFS ~ D_LAB_serum_beta2_microglobulin + HB_m_105 + 
               D_LAB_chem_ldh + MMrisk3 + sctflag, data = import)
mv6 %>% summary()

ggforest2(mv6) 
ggsave(paste0(outpath,"FOREST_plot_mmrisk3_PFS.png"), main = "PFS 1q&13-", width = 8, height = 8)




#MV7 OS UltraLow
mv7 <- coxph(OS ~ OLD + female_gender + HB_m_105 + PLT_m_150 +
               D_LAB_chem_ldh + ULTRA_Low + sctflag + `DEL_maj_focal_TP53` + SeqWGS_MAFB_CALL , data = import  )
mv7 %>% summary()

mv7 <- coxph(OS ~ female_gender + PLT_m_150 + D_LAB_chem_ldh + ULTRA_Low + sctflag , data = import  )
mv7 %>% summary()

ggforest2(mv7) 
ggsave(paste0(outpath,"FOREST_plot_UltraLow_OS.png"), main = "OS Ultra Low Risk", width = 8, height = 8)



#MV8 PFS UltraLow
mv8 <- with(import, coxph(PFS ~ OLD + female_gender + HB_m_105 + PLT_m_150 +
                            D_LAB_chem_ldh + ULTRA_Low + sctflag, data = import  )) 
mv8 %>% summary()
mv8 <- coxph(PFS ~ female_gender + HB_m_105 + PLT_m_150 +
               D_LAB_chem_ldh + ULTRA_Low + sctflag, data = import  )
mv8 %>% summary()

ggforest2(mv8) 
ggsave(paste0(outpath,"FOREST_plot_UltraLow_PFS.png"), main = "PFS Ultra Low Risk", width = 8, height = 8)


#_____ multivariate report ______
write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_190321.txt"), append = T)
write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_190321.txt"), append = T)
write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_190321.txt"), append = T)
write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_190321.txt"), append = T)
write_tsv(data.frame("Multivariate"),paste0(outpath,"report_cox_analysis_190321.txt"), append = T)

write_tsv(data.frame("MULTIVARIATA 1 - MMrisk1 - OS "),paste0(outpath,"report_cox_analysis_190321.txt"), append = T)
with(import, coxph(OS ~ OLD + female_gender + D_LAB_serum_beta2_microglobulin + D_LAB_chem_albumin + HB_m_105 + PLT_m_150 + D_LAB_chem_ldh + MMrisk1 + sctflag + `DEL_maj_focal_TP53` + SeqWGS_MAFB_CALL , data = import )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_190321.txt"), append = T, col_names = T)
write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_190321.txt"), append = T)
with(import, coxph(OS ~ female_gender + D_LAB_serum_beta2_microglobulin + PLT_m_150 + D_LAB_chem_ldh + MMrisk1 + sctflag, data = import))  %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_190321.txt"), append = T, col_names = T)



write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_190321.txt"), append = T)
write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_190321.txt"), append = T)
write_tsv(data.frame("MULTIVARIATA 2 - MMrisk1 - PFS "),paste0(outpath,"report_cox_analysis_190321.txt"), append = T)
with(import, coxph(PFS ~ OLD + female_gender + D_LAB_serum_beta2_microglobulin + D_LAB_chem_albumin + HB_m_105 + PLT_m_150 + D_LAB_chem_ldh + MMrisk1 + sctflag, data = import )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_190321.txt"), append = T, col_names = T)
write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_190321.txt"), append = T)
with(import, coxph(PFS ~  female_gender + D_LAB_serum_beta2_microglobulin + HB_m_105 + D_LAB_chem_ldh + MMrisk1 + sctflag, data = import )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_190321.txt"), append = T, col_names = T)



write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_190321.txt"), append = T)
write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_190321.txt"), append = T)
write_tsv(data.frame("MULTIVARIATA 3 - Ultra High - OS "),paste0(outpath,"report_cox_analysis_190321.txt"), append = T)
with(import, coxph(OS ~ OLD + female_gender + HB_m_105 + PLT_m_150 + D_LAB_chem_ldh + ULTRA_High + sctflag + `DEL_maj_focal_TP53` + SeqWGS_MAFB_CALL , data = import  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_190321.txt"), append = T, col_names = T)
write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_190321.txt"), append = T)
with(import, coxph(OS ~ female_gender + PLT_m_150 + D_LAB_chem_ldh + ULTRA_High + sctflag, data = import)) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_190321.txt"), append = T, col_names = T)



write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_190321.txt"), append = T)
write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_190321.txt"), append = T)
write_tsv(data.frame("MULTIVARIATA 4 - Ultra High - PFS "),paste0(outpath,"report_cox_analysis_190321.txt"), append = T)
with(import, coxph(PFS ~ OLD + female_gender + HB_m_105 + PLT_m_150 + D_LAB_chem_ldh + ULTRA_High + sctflag, data = import  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_190321.txt"), append = T, col_names = T)
write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_190321.txt"), append = T)
with(import, coxph(PFS ~ OLD + female_gender + HB_m_105 + PLT_m_150 + D_LAB_chem_ldh + ULTRA_High + sctflag, data = import )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_190321.txt"), append = T, col_names = T)



write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_190321.txt"), append = T)
write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_190321.txt"), append = T)
write_tsv(data.frame("MULTIVARIATA 5 - MMrisk3 - OS "),paste0(outpath,"report_cox_analysis_190321.txt"), append = T)
with(import, coxph(OS ~ OLD + female_gender + D_LAB_serum_beta2_microglobulin + D_LAB_chem_albumin + HB_m_105 + PLT_m_150 + D_LAB_chem_ldh + MMrisk3 + sctflag + `DEL_maj_focal_TP53` + SeqWGS_MAFB_CALL , data = import  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_190321.txt"), append = T, col_names = T)
write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_190321.txt"), append = T)
with(import, coxph(OS ~ female_gender + D_LAB_serum_beta2_microglobulin + PLT_m_150 + D_LAB_chem_ldh + MMrisk3 + sctflag , data = import )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_190321.txt"), append = T, col_names = T)



write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_190321.txt"), append = T)
write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_190321.txt"), append = T)
write_tsv(data.frame("MULTIVARIATA 6 - MMrisk3 - PFS "),paste0(outpath,"report_cox_analysis_190321.txt"), append = T)
with(import, coxph(PFS ~ OLD + female_gender + D_LAB_serum_beta2_microglobulin + D_LAB_chem_albumin + HB_m_105 + PLT_m_150 + D_LAB_chem_ldh + MMrisk3 + sctflag, data = import )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_190321.txt"), append = T, col_names = T)
write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_190321.txt"), append = T)
with(import, coxph(PFS ~ D_LAB_serum_beta2_microglobulin + HB_m_105 + D_LAB_chem_ldh + MMrisk3 + sctflag, data = import)) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_190321.txt"), append = T, col_names = T)



write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_190321.txt"), append = T)
write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_190321.txt"), append = T)
write_tsv(data.frame("MULTIVARIATA 7 - UltraLow - OS "),paste0(outpath,"report_cox_analysis_190321.txt"), append = T)
with(import, coxph(OS ~ OLD + female_gender + HB_m_105 + PLT_m_150 + D_LAB_chem_ldh + ULTRA_Low + sctflag + `DEL_maj_focal_TP53` + SeqWGS_MAFB_CALL , data = import  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_190321.txt"), append = T, col_names = T)
write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_190321.txt"), append = T)
with(import, coxph(OS ~ female_gender + PLT_m_150 + D_LAB_chem_ldh + ULTRA_Low + sctflag , data = import  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .) %>% write_tsv(paste0(outpath,"report_cox_analysis_190321.txt"), append = T, col_names = T)



write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_190321.txt"), append = T)
write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_190321.txt"), append = T)
write_tsv(data.frame("MULTIVARIATA 8 - UltraLow - PFS "),paste0(outpath,"report_cox_analysis_190321.txt"), append = T)
with(import, coxph(PFS ~ OLD + female_gender + HB_m_105 + PLT_m_150 + D_LAB_chem_ldh + ULTRA_Low + sctflag, data = import  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_190321.txt"), append = T, col_names = T)
write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_190321.txt"), append = T)
with(import, coxph(PFS ~ female_gender + HB_m_105 + PLT_m_150 + D_LAB_chem_ldh + ULTRA_Low + sctflag, data = import  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_190321.txt"), append = T, col_names = T)

