library(tidyverse)
library(survival)
library(survminer)
library(RODBC)
library(broom)


import <- data.table::fread("results/complete_database_CoMMpass_1q_13.txt") 






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

outpath <- "results/Clinical_anlysis/CoMM_dataset/"
dir.create(outpath,recursive = T, showWarnings = F)


write_tsv(data.frame("Univariate"),paste0(outpath,"report_cox_analysis_CoMMpass.txt"), append = T)



#_______ OLD ________ 

gg <- ggsurvplot(survfit(OS ~ import$OLD, data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("AGE <65 YEARS", "AGE >= 65 YEARS"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = survminer:::.build_ggsurvplot(gg), filename = "OLD_OS.png", path = outpath, dpi = 300, height = 6, width = 8)
gg <- ggsurvplot(survfit(PFS ~ import$OLD, data = import) , pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("AGE <65 YEARS", "AGE >= 65 YEARS"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = survminer:::.build_ggsurvplot(gg), filename = "OLD_PFS.png", path = outpath, dpi = 300, height = 6, width = 8)

write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_CoMMpass.txt"), append = T)
with(import, coxph(OS ~ OLD  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_CoMMpass.txt"), append = T, col_names = T)
with(import, coxph(PFS ~ OLD  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_CoMMpass.txt"), append = T, col_names = T)


#_______ SEX ________ 

gg <- ggsurvplot(survfit(OS ~ import$SEX, data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("SEX=F", "SEX=M"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = survminer:::.build_ggsurvplot(gg), filename = "SEX_OS.png", path = outpath, dpi = 300, height = 6, width = 8)
gg <- ggsurvplot(survfit(PFS ~ import$SEX, data = import) , pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("SEX=F", "SEX=M"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = survminer:::.build_ggsurvplot(gg), filename = "SEX_PFS.png", path = outpath, dpi = 300, height = 6, width = 8)

write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_CoMMpass.txt"), append = T)
with(import, coxph(OS ~ SEX )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_CoMMpass.txt"), append = T, col_names = T)
with(import, coxph(PFS ~ SEX )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_CoMMpass.txt"), append = T, col_names = T)


#_______ B2M _______ 

gg <- ggsurvplot(survfit(OS ~ import$B2M_m_median, data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("B2M high", "B2M low"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = survminer:::.build_ggsurvplot(gg), filename = "B2M_median_OS.png", path = outpath, dpi = 300, height = 6, width = 8)
gg <- ggsurvplot(survfit(PFS ~ import$HB_m_105, data = import) , pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("B2M high", "B2M low"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = survminer:::.build_ggsurvplot(gg), filename = "B2M_median_PFS.png", path = outpath, dpi = 300, height = 6, width = 8)


write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_CoMMpass.txt"), append = T)
with(import, coxph(OS ~ D_LAB_serum_beta2_microglobulin  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_CoMMpass.txt"), append = T, col_names = T)
with(import, coxph(PFS ~ D_LAB_serum_beta2_microglobulin  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_CoMMpass.txt"), append = T, col_names = T)



#_______ ALBUMIN _______ 

gg <- ggsurvplot(survfit(OS ~ import$ALBUMIN_m_median, data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("ALBUMIN high", "ALBUMIN low"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = survminer:::.build_ggsurvplot(gg), filename = "ALBUMIN_median_OS.png", path = outpath, dpi = 300, height = 6, width = 8)
gg <- ggsurvplot(survfit(PFS ~ import$ALBUMIN_m_median, data = import) , pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("ALBUMIN high", "ALBUMIN low"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = survminer:::.build_ggsurvplot(gg), filename = "ALBUMIN_median_PFS.png", path = outpath, dpi = 300, height = 6, width = 8)


write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_CoMMpass.txt"), append = T)
with(import, coxph(OS ~ D_LAB_chem_albumin  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_CoMMpass.txt"), append = T, col_names = T)
with(import, coxph(PFS ~ D_LAB_chem_albumin )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_CoMMpass.txt"), append = T, col_names = T)




#_______ HB _______ 

gg <- ggsurvplot(survfit(OS ~ import$HB_m_105, data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("HB>10.5 g/dL", "HB<10.5 g/dL"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = survminer:::.build_ggsurvplot(gg), filename = "HB_m_105_OS.png", path = outpath, dpi = 300, height = 6, width = 8)
gg <- ggsurvplot(survfit(PFS ~ import$HB_m_105, data = import) , pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("HB>10.5 g/dL", "HB<10.5 g/dL"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = survminer:::.build_ggsurvplot(gg), filename = "HB_m_105_PFS.png", path = outpath, dpi = 300, height = 6, width = 8)


write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_CoMMpass.txt"), append = T)
with(import, coxph(OS ~ HB_m_105  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_CoMMpass.txt"), append = T, col_names = T)
with(import, coxph(PFS ~ HB_m_105  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_CoMMpass.txt"), append = T, col_names = T)


#_______ PLT ________ 

gg <- ggsurvplot(survfit(OS ~ import$PLT_m_150, data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("PLT>150.000/mm³", "PLT<150.000/mm³"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = survminer:::.build_ggsurvplot(gg), filename = "PLT_m_150_OS.png", path = outpath, dpi = 300, height = 6, width = 8)
gg <- ggsurvplot(survfit(PFS ~ import$PLT_m_150, data = import) , pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("PLT>150.000/mm³", "PLT<150.000/mm³"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = survminer:::.build_ggsurvplot(gg), filename = "PLT_m_150_PFS.png", path = outpath, dpi = 300, height = 6, width = 8)


write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_CoMMpass.txt"), append = T)
with(import, coxph(OS ~ PLT_m_150  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_CoMMpass.txt"), append = T, col_names = T)
with(import, coxph(PFS ~ PLT_m_150  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_CoMMpass.txt"), append = T, col_names = T)

#________ PC ________ 

gg <- ggsurvplot(survfit(OS ~ import$PC_M_60, data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("PC<60%", "PC>60%"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = survminer:::.build_ggsurvplot(gg), filename = "PC_M_60_OS.png", path = outpath, dpi = 300, height = 6, width = 8)
gg <- ggsurvplot(survfit(PFS ~ import$PC_M_60, data = import) , pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("PC<60%", "PC>60%"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = survminer:::.build_ggsurvplot(gg), filename = "PC_M_60_PFS.png", path = outpath, dpi = 300, height = 6, width = 8)

write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_CoMMpass.txt"), append = T)
with(import, coxph(OS ~ PC_M_60  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_CoMMpass.txt"), append = T, col_names = T)
with(import, coxph(PFS ~ PC_M_60  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_CoMMpass.txt"), append = T, col_names = T)

#______ LDH __________ 

gg <- ggsurvplot(survfit(OS ~ import$LDH_level, data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = survminer:::.build_ggsurvplot(gg), filename = "LDH_UL_OS.png", path = outpath, dpi = 300, height = 6, width = 8)
gg <- ggsurvplot(survfit(PFS ~ import$LDH_level, data = import) , pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = survminer:::.build_ggsurvplot(gg), filename = "LDH_UL_PFS.png", path = outpath, dpi = 300, height = 6, width = 8)

write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_CoMMpass.txt"), append = T)
with(import, coxph(OS ~ D_LAB_chem_ldh  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_CoMMpass.txt"), append = T, col_names = T)
with(import, coxph(PFS ~ D_LAB_chem_ldh  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_CoMMpass.txt"), append = T, col_names = T)


#______ R ISS __________ 

gg <- ggsurvplot(survfit(OS ~ import$R_ISS, data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("ISS 1", "ISS 2", "ISS 3"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = survminer:::.build_ggsurvplot(gg), filename = "ISS_OS.png", path = outpath, dpi = 300, height = 6, width = 8)
gg <- ggsurvplot(survfit(PFS ~ import$R_ISS, data = import) , pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("ISS 1", "ISS 2", "ISS 3"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = survminer:::.build_ggsurvplot(gg), filename = "ISS_PFS.png", path = outpath, dpi = 300, height = 6, width = 8)

write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_CoMMpass.txt"), append = T)
with(import, coxph(OS ~ R_ISS  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_CoMMpass.txt"), append = T, col_names = T)
with(import, coxph(PFS ~ R_ISS  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_CoMMpass.txt"), append = T, col_names = T)

#______ Light-Chain __________ 

gg <- ggsurvplot(survfit(OS ~ import$D_IM_LIGHT_CHAIN_BY_FLOW, data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="",  tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = survminer:::.build_ggsurvplot(gg), filename = "LIGHT_CHAIN_OS.png", path = outpath, dpi = 300, height = 6, width = 8)
gg <- ggsurvplot(survfit(PFS ~ import$D_IM_LIGHT_CHAIN_BY_FLOW, data = import) , pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="",  tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = survminer:::.build_ggsurvplot(gg), filename = "LIGHT_CHAIN_PFS.png", path = outpath, dpi = 300, height = 6, width = 8)


write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_CoMMpass.txt"), append = T)
with(import, coxph(OS ~ D_IM_LIGHT_CHAIN_BY_FLOW  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_CoMMpass.txt"), append = T, col_names = T)
with(import, coxph(PFS ~ D_IM_LIGHT_CHAIN_BY_FLOW  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_CoMMpass.txt"), append = T, col_names = T)


#______ TX yes/no ________ 

gg <- ggsurvplot(survfit(OS ~ import$sctflag, data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("TX_I_LINE= 0", "TX_I_LINE= 1"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = survminer:::.build_ggsurvplot(gg), filename = "TX_I_LINE_1_0_OS.png", path = outpath, dpi = 300, height = 6, width = 8)
gg <- ggsurvplot(survfit(PFS ~ import$sctflag, data = import) , pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("TX_I_LINE= 0", "TX_I_LINE= 1"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = survminer:::.build_ggsurvplot(gg), filename = "TX_I_LINE_1_0_PFS.png", path = outpath, dpi = 300, height = 6, width = 8)

write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_CoMMpass.txt"), append = T)
with(import, coxph(OS ~ sctflag  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_CoMMpass.txt"), append = T, col_names = T)
with(import, coxph(PFS ~ sctflag  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_CoMMpass.txt"), append = T, col_names = T)

#_______ t(4;14) _______ 

gg <- ggsurvplot(survfit(OS ~ import$SeqWGS_WHSC1_CALL, data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("FISH_T_4_14=0", "FISH_T_4_14=1" ), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = survminer:::.build_ggsurvplot(gg), filename = "FISH_T_4_14_OS.png", path = outpath, dpi = 300, height = 6, width = 8)
gg <- ggsurvplot(survfit(PFS ~ import$SeqWGS_WHSC1_CALL, data = import) , pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("FISH_T_4_14=0", "FISH_T_4_14=1"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = survminer:::.build_ggsurvplot(gg), filename = "FISH_T_4_14_PFS.png", path = outpath, dpi = 300, height = 6, width = 8)

write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_CoMMpass.txt"), append = T)
with(import, coxph(OS ~ SeqWGS_WHSC1_CALL  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_CoMMpass.txt"), append = T, col_names = T)
with(import, coxph(PFS ~ SeqWGS_WHSC1_CALL  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_CoMMpass.txt"), append = T, col_names = T)

#_______ t(11;14) _______ 

gg <- ggsurvplot(survfit(OS ~ import$SeqWGS_CCND1_CALL, data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("FISH_T_11_14=0", "FISH_T_11_14=1" ), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = survminer:::.build_ggsurvplot(gg), filename = "FISH_T_11_14_OS.png", path = outpath, dpi = 300, height = 6, width = 8)
gg <- ggsurvplot(survfit(PFS ~ import$SeqWGS_CCND1_CALL, data = import) , pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("FISH_T_11_14=0", "FISH_T_11_14=1"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = survminer:::.build_ggsurvplot(gg), filename = "FISH_T_11_14_PFS.png", path = outpath, dpi = 300, height = 6, width = 8)

write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_CoMMpass.txt"), append = T)
with(import, coxph(OS ~ SeqWGS_CCND1_CALL  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_CoMMpass.txt"), append = T, col_names = T)
with(import, coxph(PFS ~ SeqWGS_CCND1_CALL  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_CoMMpass.txt"), append = T, col_names = T)

#_______ t(14;16) _______ 

gg <- ggsurvplot(survfit(OS ~ import$SeqWGS_MAF_CALL, data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("FISH_T_14_16=0", "FISH_T_14_16=1" ), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = survminer:::.build_ggsurvplot(gg), filename = "FISH_T_14_16_OS.png", path = outpath, dpi = 300, height = 6, width = 8)
gg <- ggsurvplot(survfit(PFS ~ import$SeqWGS_MAF_CALL, data = import) , pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("FISH_T_14_16=0", "FISH_T_14_16=1"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = survminer:::.build_ggsurvplot(gg), filename = "FISH_T_14_16_PFS.png", path = outpath, dpi = 300, height = 6, width = 8)

write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_CoMMpass.txt"), append = T)
with(import, coxph(OS ~ SeqWGS_MAF_CALL  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_CoMMpass.txt"), append = T, col_names = T)
with(import, coxph(PFS ~ SeqWGS_MAF_CALL  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_CoMMpass.txt"), append = T, col_names = T)

#_______ t(14;20) _______ 

gg <- ggsurvplot(survfit(OS ~ import$SeqWGS_MAFB_CALL, data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("FISH_T_14_20=0", "FISH_T_14_20=1" ), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = survminer:::.build_ggsurvplot(gg), filename = "FISH_T_14_20_OS.png", path = outpath, dpi = 300, height = 6, width = 8)
gg <- ggsurvplot(survfit(PFS ~ import$SeqWGS_MAFB_CALL, data = import) , pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("FISH_T_14_20=0", "FISH_T_14_20=1"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = survminer:::.build_ggsurvplot(gg), filename = "FISH_T_14_20_PFS.png", path = outpath, dpi = 300, height = 6, width = 8)

write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_CoMMpass.txt"), append = T)
with(import, coxph(OS ~ SeqWGS_MAFB_CALL  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_CoMMpass.txt"), append = T, col_names = T)
with(import, coxph(PFS ~ SeqWGS_MAFB_CALL  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_CoMMpass.txt"), append = T, col_names = T)

#_______ t(6;14) _______ 

gg <- ggsurvplot(survfit(OS ~ import$SeqWGS_CCND3_CALL, data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("FISH_T_6_14=0", "FISH_T_6_14=1" ), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = survminer:::.build_ggsurvplot(gg), filename = "FISH_T_6_14_OS.png", path = outpath, dpi = 300, height = 6, width = 8)
gg <- ggsurvplot(survfit(PFS ~ import$SeqWGS_CCND3_CALL, data = import) , pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("FISH_T_6_14=0", "FISH_T_6_14=1"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = survminer:::.build_ggsurvplot(gg), filename = "FISH_T_6_14_PFS.png", path = outpath, dpi = 300, height = 6, width = 8)

write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_CoMMpass.txt"), append = T)
with(import, coxph(OS ~ SeqWGS_CCND3_CALL  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_CoMMpass.txt"), append = T, col_names = T)
with(import, coxph(PFS ~ SeqWGS_CCND3_CALL  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_CoMMpass.txt"), append = T, col_names = T)

#_______ Hyperdiploidy _______ 

gg <- ggsurvplot(survfit(OS ~ import$Hyperdiploid, data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("HyperDiploidy=0", "HyperDiploidy=1" ), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = survminer:::.build_ggsurvplot(gg), filename = "HyperDiploidy_OS.png", path = outpath, dpi = 300, height = 6, width = 8)
gg <- ggsurvplot(survfit(PFS ~ import$Hyperdiploid, data = import) , pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("HyperDiploidy=0", "HyperDiploidy=1"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = survminer:::.build_ggsurvplot(gg), filename = "HyperDiploidy_PFS.png", path = outpath, dpi = 300, height = 6, width = 8)


write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_CoMMpass.txt"), append = T)
with(import, coxph(OS ~ Hyperdiploid  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_CoMMpass.txt"), append = T, col_names = T)
with(import, coxph(PFS ~ Hyperdiploid  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_CoMMpass.txt"), append = T, col_names = T)


#______ DEL 17p _______ 

gg <- ggsurvplot(survfit(OS ~ import$"DEL_maj_broad_chr_17p", data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("Chr 17p Del= 0", "Chr 17p Del= 1" ), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = survminer:::.build_ggsurvplot(gg), filename = "DEL_maj_broad_chr_17p_OS.png", path = outpath, dpi = 300, height = 6, width = 8)
gg <- ggsurvplot(survfit(PFS ~ import$"DEL_maj_broad_chr_17p", data = import) , pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("Chr 17p Del= 0", "Chr 17p Del= 1"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = survminer:::.build_ggsurvplot(gg), filename = "DEL_maj_broad_chr_17p_PFS.png", path = outpath, dpi = 300, height = 6, width = 8)

write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_CoMMpass.txt"), append = T)
with(import, coxph(OS ~ `DEL_maj_broad_chr_17p`  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_CoMMpass.txt"), append = T, col_names = T)
with(import, coxph(PFS ~ `DEL_maj_broad_chr_17p`  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_CoMMpass.txt"), append = T, col_names = T)


#______ DEL 17p _______ 

gg <- ggsurvplot(survfit(OS ~ import$"DEL_maj_broad_chr_1p", data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("Chr 1p Del= 0", "Chr 1p Del= 1" ), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = survminer:::.build_ggsurvplot(gg), filename = "DEL_maj_broad_chr_1p_OS.png", path = outpath, dpi = 300, height = 6, width = 8)
gg <- ggsurvplot(survfit(PFS ~ import$"DEL_maj_broad_chr_1p", data = import) , pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("Chr 1p Del= 0", "Chr 1p Del= 1"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = survminer:::.build_ggsurvplot(gg), filename = "DEL_maj_broad_chr_1p_PFS.png", path = outpath, dpi = 300, height = 6, width = 8)

write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_CoMMpass.txt"), append = T)
with(import, coxph(OS ~ `DEL_maj_broad_chr_1p`  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_CoMMpass.txt"), append = T, col_names = T)
with(import, coxph(PFS ~ `DEL_maj_broad_chr_1p`  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_CoMMpass.txt"), append = T, col_names = T)



#______ DEL TP53 ________ 

gg <- ggsurvplot(survfit(OS ~ import$"DEL_maj_focal_TP53", data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("TP53 Del= 0", "TP53 Del= 1" ), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = survminer:::.build_ggsurvplot(gg), filename = "DEL_maj_focal_TP53_OS.png", path = outpath, dpi = 300, height = 6, width = 8)
gg <- ggsurvplot(survfit(PFS ~ import$"DEL_maj_focal_TP53", data = import) , pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("TP53 Del= 0", "TP53 Del= 1"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = survminer:::.build_ggsurvplot(gg), filename = "DEL_maj_focal_TP53_PFS.png", path = outpath, dpi = 300, height = 6, width = 8)

write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_CoMMpass.txt"), append = T)
with(import, coxph(OS ~ `DEL_maj_focal_TP53`  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_CoMMpass.txt"), append = T, col_names = T)
with(import, coxph(PFS ~ `DEL_maj_focal_TP53`  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_CoMMpass.txt"), append = T, col_names = T)

#______ AMP 1q ________ 

gg <- ggsurvplot(survfit(OS ~ import$"AMP_maj_broad_chr_1q", data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("Chr 1q Broad Amp= 0", "Chr 1q Broad Amp= 1" ), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = survminer:::.build_ggsurvplot(gg), filename = "AMP_maj_broad_chr_1_OS.png", path = outpath, dpi = 300, height = 6, width = 8)
gg <- ggsurvplot(survfit(PFS ~ import$"AMP_maj_broad_chr_1q", data = import) , pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("Chr 1q Broad Amp= 0", "Chr 1q Broad Amp= 1"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = survminer:::.build_ggsurvplot(gg), filename = "AMP_maj_broad_chr_1_PFS.png", path = outpath, dpi = 300, height = 6, width = 8)

write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_CoMMpass.txt"), append = T)
with(import, coxph(OS ~ `AMP_maj_broad_chr_1q`  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_CoMMpass.txt"), append = T, col_names = T)
with(import, coxph(PFS ~ `AMP_maj_broad_chr_1q`  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_CoMMpass.txt"), append = T, col_names = T)



# 1q MMrisk focal+broad 

gg <- ggsurvplot(survfit(OS ~ import$MMrisk_1q_all, data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("Chr1q Broad&Focal Amp=0", "Chr1q Broad&Focal Amp=1" ), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = survminer:::.build_ggsurvplot(gg), filename = "MMrisk_1q_all_OS.png", path = outpath, dpi = 300, height = 6, width = 8)
gg <- ggsurvplot(survfit(PFS ~ import$MMrisk_1q_all, data = import) , pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("Chr1q Broad&Focal Amp=0", "Chr1q Broad&Focal Amp=1"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = survminer:::.build_ggsurvplot(gg), filename = "MMrisk_1q_all_PFS.png", path = outpath, dpi = 300, height = 6, width = 8)

write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_CoMMpass.txt"), append = T)
with(import, coxph(OS ~ MMrisk_1q_all  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_CoMMpass.txt"), append = T, col_names = T)
with(import, coxph(PFS ~ MMrisk_1q_all  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_CoMMpass.txt"), append = T, col_names = T)


#______ DEL 13 ________ 

gg <- ggsurvplot(survfit(OS ~ import$"DEL_maj_broad_chr_13q", data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("Chr 13q Broad Del= 0", "Chr 13q Broad Del= 1" ), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = survminer:::.build_ggsurvplot(gg), filename = "DEL_maj_broad_chr_13q_OS.png", path = outpath, dpi = 300, height = 6, width = 8)
gg <- ggsurvplot(survfit(PFS ~ import$"DEL_maj_broad_chr_13q", data = import) , pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("Chr 13q Broad Del= 0", "Chr 13q Broad Del= 1"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = survminer:::.build_ggsurvplot(gg), filename = "DEL_maj_broad_chr_13q_PFS.png", path = outpath, dpi = 300, height = 6, width = 8)


write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_CoMMpass.txt"), append = T)
with(import, coxph(OS ~ `DEL_maj_broad_chr_13q`  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_CoMMpass.txt"), append = T, col_names = T)
with(import, coxph(PFS ~ `DEL_maj_broad_chr_13q`  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_CoMMpass.txt"), append = T, col_names = T)


# 13 MMrisk focal+broad

gg <- ggsurvplot(survfit(OS ~ import$MMrisk_13_all, data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("Chr 13q Broad&Focal Del= 0", "Chr 13q Broad&Focal Del= 1" ), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = survminer:::.build_ggsurvplot(gg), filename = "MMrisk_13_all_OS.png", path = outpath, dpi = 300, height = 6, width = 8)
gg <- ggsurvplot(survfit(PFS ~ import$MMrisk_13_all, data = import) , pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("Chr 13q Broad&Focal Del= 0", "Chr 13q Broad&Focal Del= 1"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = survminer:::.build_ggsurvplot(gg), filename = "MMrisk_13_all_PFS.png", path = outpath, dpi = 300, height = 6, width = 8)

write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_CoMMpass.txt"), append = T)
with(import, coxph(OS ~ MMrisk_13_all  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_CoMMpass.txt"), append = T, col_names = T)
with(import, coxph(PFS ~ MMrisk_13_all  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_CoMMpass.txt"), append = T, col_names = T)

#______ MM RISK ________ 

gg <- ggsurvplot(survfit(OS ~ import$MMrisk_CLASSIFIER_ALL, data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("1q&13+", "1q/13", "1q&13-" ), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = survminer:::.build_ggsurvplot(gg), filename = "MMrisk_CLASSIFIER_ALL_OS.png", path = outpath, dpi = 300, height = 6, width = 8)
gg <- ggsurvplot(survfit(PFS ~ import$MMrisk_CLASSIFIER_ALL, data = import) , pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("1q&13+", "1q/13","1q&13-" ), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = survminer:::.build_ggsurvplot(gg), filename = "MMrisk_CLASSIFIER_ALL_PFS.png", path = outpath, dpi = 300, height = 6, width = 8)


write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_CoMMpass.txt"), append = T)
with(import, coxph(OS ~ MMrisk_CLASSIFIER_ALL  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_CoMMpass.txt"), append = T, col_names = T)
with(import, coxph(PFS ~ MMrisk_CLASSIFIER_ALL  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_CoMMpass.txt"), append = T, col_names = T)

survfit(OS ~ MMrisk_CLASSIFIER_ALL, data=import  ) %>% summary
survfit(PFS ~ MMrisk_CLASSIFIER_ALL, data=import  ) %>% summary



#______ Risk 1 ______ 

gg <- ggsurvplot(survfit(OS ~ import$MMrisk1, data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("1q&13+ =0", "1q&13+ =1" ), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = survminer:::.build_ggsurvplot(gg), filename = "MMrisk1_OS.png", path = outpath, dpi = 300, height = 6, width = 8)
gg <- ggsurvplot(survfit(PFS ~ import$MMrisk1, data = import) , pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("1q&13+ =0", "1q&13+ =1"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = survminer:::.build_ggsurvplot(gg), filename = "MMrisk1_PFS.png", path = outpath, dpi = 300, height = 6, width = 8)

write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_CoMMpass.txt"), append = T)
with(import, coxph(OS ~ MMrisk1  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_CoMMpass.txt"), append = T, col_names = T)
with(import, coxph(PFS ~ MMrisk1  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_CoMMpass.txt"), append = T, col_names = T)

#_______ Risk 2 ______ 

gg <- ggsurvplot(survfit(OS ~ import$MMrisk2, data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("1q/13+ =0", "1q/13+ =1" ), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = survminer:::.build_ggsurvplot(gg), filename = "MMrisk2_OS.png", path = outpath, dpi = 300, height = 6, width = 8)
gg <- ggsurvplot(survfit(PFS ~ import$MMrisk2, data = import) , pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("1q/13+ =0", "1q/13+ =1"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = survminer:::.build_ggsurvplot(gg), filename = "MMrisk2_PFS.png", path = outpath, dpi = 300, height = 6, width = 8)


write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_CoMMpass.txt"), append = T)
with(import, coxph(OS ~ MMrisk2  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_CoMMpass.txt"), append = T, col_names = T)
with(import, coxph(PFS ~ MMrisk2  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_CoMMpass.txt"), append = T, col_names = T)


#_______ Risk 3 _______ 

gg <- ggsurvplot(survfit(OS ~ import$MMrisk3, data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("1q&13- =0", "1q&13- =1" ), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = survminer:::.build_ggsurvplot(gg), filename = "MMrisk3_OS.png", path = outpath, dpi = 300, height = 6, width = 8)
gg <- ggsurvplot(survfit(PFS ~ import$MMrisk3, data = import) , pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("1q&13- =0", "1q&13- =1"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = survminer:::.build_ggsurvplot(gg), filename = "MMrisk3_PFS.png", path = outpath, dpi = 300, height = 6, width = 8)


write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_CoMMpass.txt"), append = T)
with(import, coxph(OS ~ MMrisk3  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_CoMMpass.txt"), append = T, col_names = T)
with(import, coxph(PFS ~ MMrisk3  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_CoMMpass.txt"), append = T, col_names = T)


#____ MM RISK a/b _______

gg <- ggsurvplot(survfit(OS ~ import$MMrisk_AB_ALL, data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("1q&13+", "only 13 Del", "only 1q Amp", "1q&13-" ), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = survminer:::.build_ggsurvplot(gg), filename = "MMrisk_AB_ALL_OS.png", path = outpath, dpi = 300, height = 6, width = 8)
gg <- ggsurvplot(survfit(PFS ~ import$MMrisk_AB_ALL, data = import) , pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("1q&13+", "only 13 Del", "only 1q Amp", "1q&13-"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = survminer:::.build_ggsurvplot(gg), filename = "MMrisk_AB_ALL_PFS.png", path = outpath, dpi = 300, height = 6, width = 8)



write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_CoMMpass.txt"), append = T)
with(import, coxph(OS ~ MMrisk_AB_ALL  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_CoMMpass.txt"), append = T, col_names = T)
with(import, coxph(PFS ~ MMrisk_AB_ALL  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_CoMMpass.txt"), append = T, col_names = T)



#____ Risk_1 vs only_1q _______

gg <- ggsurvplot(survfit(OS ~ import$risk1_VS_1q_VS_other, data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("1q&13+", "only 1q Amp", "other" ), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = survminer:::.build_ggsurvplot(gg), filename = "risk1_VS_1q_VS_other_OS.png", path = outpath, dpi = 300, height = 6, width = 8)
gg <- ggsurvplot(survfit(PFS ~ import$risk1_VS_1q_VS_other, data = import) , pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("1q&13+", "only 1q amp", "other"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = survminer:::.build_ggsurvplot(gg), filename = "risk1_VS_1q_VS_other_PFS.png", path = outpath, dpi = 300, height = 6, width = 8)



write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_CoMMpass.txt"), append = T)
with(import, coxph(OS ~ risk1_VS_1q_VS_other  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_CoMMpass.txt"), append = T, col_names = T)
with(import, coxph(PFS ~ risk1_VS_1q_VS_other  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_CoMMpass.txt"), append = T, col_names = T)



imp1q <- filter(import, import$"AMP_maj_broad_chr_1q" == 1)

OS1q <- Surv(imp1q$OS_months, imp1q$censos)
PFS1q <- Surv(imp1q$PFS_months, imp1q$censpfs)

gg <- ggsurvplot(survfit(OS1q ~ imp1q$risk1_VS_1q_VS_other, data = imp1q) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("1q&13+", "only 1q Amp" ), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = survminer:::.build_ggsurvplot(gg), filename = "risk1_VS_1q_VS_other_OS1q.png", path = outpath, dpi = 300, height = 6, width = 8)
gg <- ggsurvplot(survfit(PFS1q ~ imp1q$risk1_VS_1q_VS_other, data = imp1q) , pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("1q&13+", "only 1q Amp"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = survminer:::.build_ggsurvplot(gg), filename = "risk1_VS_1q_VS_other_PFS1q.png", path = outpath, dpi = 300, height = 6, width = 8)


write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_CoMMpass.txt"), append = T)
with(imp1q, coxph(OS1q ~ risk1_VS_1q_VS_other  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_CoMMpass.txt"), append = T, col_names = T)
with(imp1q, coxph(PFS1q ~ risk1_VS_1q_VS_other  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_CoMMpass.txt"), append = T, col_names = T)


#______ ultra MMrisk classificator ______

gg <- ggsurvplot(survfit(OS ~ import$ULTRA_MMrisk, data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("ULTRA_MMrisk=1", "ULTRA_MMrisk=2", "ULTRA_MMrisk=3" ), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = survminer:::.build_ggsurvplot(gg), filename = "ULTRA_MMrisk_OS.png", path = outpath, dpi = 300, height = 6, width = 8)
gg <- ggsurvplot(survfit(PFS ~ import$ULTRA_MMrisk, data = import), pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("ULTRA_MMrisk=1", "ULTRA_MMrisk=2", "ULTRA_MMrisk=3"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = survminer:::.build_ggsurvplot(gg), filename = "ULTRA_MMrisk_PFS.png", path = outpath, dpi = 300, height = 6, width = 8)

write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_CoMMpass.txt"), append = T)
with(import, coxph(OS ~ ULTRA_MMrisk  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_CoMMpass.txt"), append = T, col_names = T)
with(import, coxph(PFS ~ ULTRA_MMrisk  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_CoMMpass.txt"), append = T, col_names = T)



#______ ultra high  ______

gg <- ggsurvplot(survfit(OS ~ import$ULTRA_High, data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("ULTRA_High=0", "ULTRA_High=1"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = survminer:::.build_ggsurvplot(gg), filename = "ULTRA_High_OS.png", path = outpath, dpi = 300, height = 6, width = 8)
gg <- ggsurvplot(survfit(PFS ~ import$ULTRA_High, data = import), pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("ULTRA_High=0", "ULTRA_High=1"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = survminer:::.build_ggsurvplot(gg), filename = "ULTRA_High_PFS.png", path = outpath, dpi = 300, height = 6, width = 8)

write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_CoMMpass.txt"), append = T)
with(import, coxph(OS ~ ULTRA_High  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_CoMMpass.txt"), append = T, col_names = T)
with(import, coxph(PFS ~ ULTRA_High  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_CoMMpass.txt"), append = T, col_names = T)



#______ ultra low  ______

gg <- ggsurvplot(survfit(OS ~ import$ULTRA_Low, data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("ULTRA_Low=0", "ULTRA_Low=1"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = survminer:::.build_ggsurvplot(gg), filename = "ULTRA_High_OS.png", path = outpath, dpi = 300, height = 6, width = 8)
gg <- ggsurvplot(survfit(PFS ~ import$ULTRA_Low, data = import), pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("ULTRA_Low=0", "ULTRA_Low=1"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = survminer:::.build_ggsurvplot(gg), filename = "ULTRA_High_PFS.png", path = outpath, dpi = 300, height = 6, width = 8)

write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_CoMMpass.txt"), append = T)
with(import, coxph(OS ~ ULTRA_Low  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_CoMMpass.txt"), append = T, col_names = T)
with(import, coxph(PFS ~ ULTRA_Low  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_CoMMpass.txt"), append = T, col_names = T)





################################# KM curves with ISS and R ISS #################################

import$MMrisk_class <- import$MMrisk_class %>% recode_factor(risk_1 ="1q&13+", risk_2 ="1q/13", risk_3 ="1q&13-") %>% relevel("1q&13-")
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


import$CCND2_traslocation <- 0
import$CCND2_traslocation[import$SeqFISH_T_4_14==1] <- 1
import$CCND2_traslocation[import$SeqFISH_T_14_16==1] <- 1
import$CCND2_traslocation[import$SeqFISH_T_14_20==1] <- 1

import$MMrisk_class_t_CCND2 <- ifelse(import$MMrisk_CLASS==1 & import$CCND2_traslocation ==1 , "t&1q&13+",
                                      ifelse(import$MMrisk_CLASS==1, "1q&13+_pure", import$MMrisk_class %>% as.character))

import$MMrisk_class_t_CCND2 %>% table



gg <- ggsurvplot(survfit(OS ~ import$MMrisk_class, data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)


import$MM_ISS <- paste0(import$MMrisk_class,"/R-ISS_",import$R_ISS)


gg <- ggsurvplot(survfit(OS ~ import$MM_ISS, data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)


import$MMT_ISS <- paste0(import$MMrisk_class_t_CCND2,"/R-ISS_",import$R_ISS)
import$MMT_ISS %>% table

gg <- ggsurvplot(survfit(OS ~ import$MMT_ISS, data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.5), legend.title="", tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)


i2 <- import %>% filter(R_ISS != "NA")

OS2 <- Surv(i2$OS_months, i2$censos)

gg <- ggsurvplot(survfit(OS2 ~ MMT_ISS, data = i2) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.5), legend.title="", tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)




i3 <- import %>% filter(R_ISS != "NA" & MMrisk_class_t_CCND2=="1q&13+_pure")

OS3 <- Surv(i3$OS_months, i3$censos)

gg <- ggsurvplot(survfit(OS3 ~ MMT_ISS, data = i3) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.5), legend.title="", tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)




i4 <- import %>% filter(R_ISS != "NA" & MMrisk_class_t_CCND2=="t&1q&13+")

OS4 <- Surv(i4$OS_months, i4$censos)

gg <- ggsurvplot(survfit(OS4 ~ MMT_ISS, data = i4) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.5), legend.title="", tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)



i5 <- import %>% filter(R_ISS == "3")

OS5 <- Surv(i5$OS_months, i5$censos)

gg <- ggsurvplot(survfit(OS5 ~ MMT_ISS, data = i5) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.5), legend.title="", tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)


import$FINAL <- ifelse(import$MMT_ISS %in% c("t&1q&13+/R-ISS_3", "1q&13+_pure/R-ISS_3" ), "DEF_RISK", ifelse(import$MMT_ISS %in% c("1q&13-/R-ISS_1", "1q&13+_pure/R-ISS_1"), "LOW_RISK", "other"))

import$FINAL %>% table

gg <- ggsurvplot(survfit(OS ~ FINAL, data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.5), legend.title="", tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)

gg <- ggsurvplot(survfit(PFS ~ FINAL, data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.5), legend.title="", tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)




gg <- ggsurvplot(survfit(OS ~ R_ISS, data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.5), legend.title="", tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)

gg <- ggsurvplot(survfit(PFS ~ R_ISS, data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.5), legend.title="", tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
