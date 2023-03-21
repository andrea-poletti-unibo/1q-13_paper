library(data.table)
library(tidyverse)
library(survival)
library(survminer)
library(RODBC)
library(broom)

db_clin <- RODBC::odbcConnectAccess2007("C:/Users/mm_gr/Alma Mater Studiorum Università di Bologna/MM group - Documenti/EVERGREEN_inter-projects_files/MM_db_clinical_records.accdb")
clin <- sqlFetch(db_clin, "Extraction_1q13_301120", dec=",", na.strings="nv" )

db_exp<- RODBC::odbcConnectAccess2007("C:/Users/mm_gr/Alma Mater Studiorum Università di Bologna/MM group - Documenti/EVERGREEN_inter-projects_files/MM_db_experiments.accdb")
calls <- sqlFetch(db_exp, "SNP_calls_1q_13_181019", dec=",")

import <- inner_join(clin, calls, by=c("SNP_ARRAY_ESORDIO"="SNP"))

import$PROTOCOL_REV <- recode(import$PROTOCOL, "FP-EMN02"="AMB", "FP - BO2005"="AMB", "DARA"="AMB", "FORTE"="AMB")
import$IG_ISOTYPE_REV <- recode(import$IG_ISOTYPE, "IgA+BJ"="IgA", "IgA/BJ"="IgA")


#______ MM RISK CREATION _____
import$AMP_1q_genes_all <- with(import, 
                                ifelse(`AMP_maj-focal_ANP32E`==1 | `AMP_maj-focal_MCL1` ==1 | `AMP_maj-focal_CKS1B`==1, 1,0 ))

import$MMrisk_1q_all <- ifelse( import$`AMP_maj-broad_chr_1q` == 1 | import$AMP_1q_genes_all ==1, 1, 0)
import$MMrisk_13_all <- ifelse( import$`DEL_maj-broad_chr_13q` ==1 | import$`DEL_maj-focal_RB1` ==1, 1, 0)

import$MMrisk_CLASSIFIER_ALL <- 3  - import$MMrisk_1q_all - import$MMrisk_13_all

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
import$ULTRA_MMrisk %>% table

import$ULTRA_High <- ifelse(import$ULTRA_MMrisk=="Ultra_High",1,0)
import$ULTRA_Low  <- ifelse(import$ULTRA_MMrisk=="Ultra_Low",1,0)



#______ create surv ______

OS <- Surv(import$OS_MONTHS, import$OS_EVENT)
PFS <- Surv(import$PFS_I_MONTHS, import$PFS_I_EVENT)


#______ create additional variables ______

B2M_med <- import$B2M %>% median(na.rm = T)
import$B2M_m_median <- ifelse(import$B2M>B2M_med, 1,0)

albumin_med <- import$ALBUMIN %>% median(na.rm = T)
import$ALBUMIN_m_median <- ifelse(import$ALBUMIN>albumin_med, 1,0)


################################################################################ 
################################ SURV ANALYSIS #################################
################################################################################ 


outpath <- "C:/Users/mm_gr/Alma Mater Studiorum Università di Bologna/PROJECT 1q & 13 - Documenti/Survival_analysis/analysis_240221/"


#_______ PROTOCOL ________ 

gg <- ggsurvplot(survfit(OS ~ import$PROTOCOL_REV, data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("AMB", "BO2005", "EMN02"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = print(gg), filename = "protocol_rev_OS.png", path = outpath, dpi = 300, height = 6, width = 8)
gg <- ggsurvplot(survfit(PFS ~ import$PROTOCOL_REV, data = import) , pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("AMB", "BO2005", "EMN02"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = print(gg), filename = "protocol_rev_PFS.png", path = outpath, dpi = 300, height = 6, width = 8)

with(import, coxph(OS ~ PROTOCOL_REV  )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_240221.txt"), append = T, col_names = T)
with(import, coxph(PFS ~ PROTOCOL_REV )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_240221.txt"), append = T, col_names = T)



#_______ OLD ________ 

gg <- ggsurvplot(survfit(OS ~ import$OLD_0_1, data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("AGE <65 YEARS", "AGE >= 65 YEARS"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = print(gg), filename = "OLD_0_1_OS.png", path = outpath, dpi = 300, height = 6, width = 8)
gg <- ggsurvplot(survfit(PFS ~ import$OLD_0_1, data = import) , pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("AGE <65 YEARS", "AGE >= 65 YEARS"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = print(gg), filename = "OLD_0_1_PFS.png", path = outpath, dpi = 300, height = 6, width = 8)

write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_240221.txt"), append = T)
with(import, coxph(OS ~ OLD_0_1 + strata(PROTOCOL_REV) )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_240221.txt"), append = T, col_names = T)
with(import, coxph(PFS ~ OLD_0_1 + strata(PROTOCOL_REV) )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_240221.txt"), append = T, col_names = T)


#_______ SEX ________ 

gg <- ggsurvplot(survfit(OS ~ import$SEX, data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("SEX=F", "SEX=M"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = print(gg), filename = "SEX_OS.png", path = outpath, dpi = 300, height = 6, width = 8)
gg <- ggsurvplot(survfit(PFS ~ import$SEX, data = import) , pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("SEX=F", "SEX=M"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = print(gg), filename = "SEX_PFS.png", path = outpath, dpi = 300, height = 6, width = 8)

write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_240221.txt"), append = T)
with(import, coxph(OS ~ SEX + strata(PROTOCOL_REV) )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_240221.txt"), append = T, col_names = T)
with(import, coxph(PFS ~ SEX + strata(PROTOCOL_REV) )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_240221.txt"), append = T, col_names = T)


#_______ B2M _______ 

gg <- ggsurvplot(survfit(OS ~ import$B2M_m_median, data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("B2M high", "B2M low"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = print(gg), filename = "B2M_median_OS.png", path = outpath, dpi = 300, height = 6, width = 8)
gg <- ggsurvplot(survfit(PFS ~ import$HB_m_105, data = import) , pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("B2M high", "B2M low"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = print(gg), filename = "B2M_median_PFS.png", path = outpath, dpi = 300, height = 6, width = 8)


write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_240221.txt"), append = T)
with(import, coxph(OS ~ B2M + strata(PROTOCOL_REV) )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_240221.txt"), append = T, col_names = T)
with(import, coxph(PFS ~ B2M + strata(PROTOCOL_REV) )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_240221.txt"), append = T, col_names = T)



#_______ ALBUMIN _______ 

gg <- ggsurvplot(survfit(OS ~ import$ALBUMIN_m_median, data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("ALBUMIN high", "ALBUMIN low"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = print(gg), filename = "ALBUMIN_median_OS.png", path = outpath, dpi = 300, height = 6, width = 8)
gg <- ggsurvplot(survfit(PFS ~ import$ALBUMIN_m_median, data = import) , pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("ALBUMIN high", "ALBUMIN low"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = print(gg), filename = "ALBUMIN_median_PFS.png", path = outpath, dpi = 300, height = 6, width = 8)


write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_240221.txt"), append = T)
with(import, coxph(OS ~ ALBUMIN + strata(PROTOCOL_REV) )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_240221.txt"), append = T, col_names = T)
with(import, coxph(PFS ~ ALBUMIN + strata(PROTOCOL_REV) )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_240221.txt"), append = T, col_names = T)




#_______ HB _______ 

gg <- ggsurvplot(survfit(OS ~ import$HB_m_105, data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("HB > 10.5 g/dL", "HB < 10.5 g/dL"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = print(gg), filename = "HB_m_105_OS.png", path = outpath, dpi = 300, height = 6, width = 8)
gg <- ggsurvplot(survfit(PFS ~ import$HB_m_105, data = import) , pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("HB > 10.5 g/dL", "HB < 10.5 g/dL"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = print(gg), filename = "HB_m_105_PFS.png", path = outpath, dpi = 300, height = 6, width = 8)


write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_240221.txt"), append = T)
with(import, coxph(OS ~ HB_m_105 + strata(PROTOCOL_REV) )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_240221.txt"), append = T, col_names = T)
with(import, coxph(PFS ~ HB_m_105 + strata(PROTOCOL_REV) )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_240221.txt"), append = T, col_names = T)


#_______ PLT ________ 

gg <- ggsurvplot(survfit(OS ~ import$PLT_m_150, data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("PLT> 150.000/mm³", "PLT< 150.000/mm³"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = print(gg), filename = "PLT_m_150_OS.png", path = outpath, dpi = 300, height = 6, width = 8)
gg <- ggsurvplot(survfit(PFS ~ import$PLT_m_150, data = import) , pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("PLT> 150.000/mm³", "PLT< 150.000/mm³"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = print(gg), filename = "PLT_m_150_PFS.png", path = outpath, dpi = 300, height = 6, width = 8)


write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_240221.txt"), append = T)
with(import, coxph(OS ~ PLT_m_150 + strata(PROTOCOL_REV) )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_240221.txt"), append = T, col_names = T)
with(import, coxph(PFS ~ PLT_m_150 + strata(PROTOCOL_REV) )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_240221.txt"), append = T, col_names = T)

#________ PC ________ 

gg <- ggsurvplot(survfit(OS ~ import$PC_M_60, data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("PC< 60%", "PC> 60%"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = print(gg), filename = "PC_M_60_OS.png", path = outpath, dpi = 300, height = 6, width = 8)
gg <- ggsurvplot(survfit(PFS ~ import$PC_M_60, data = import) , pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("PC< 60%", "PC> 60%"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = print(gg), filename = "PC_M_60_PFS.png", path = outpath, dpi = 300, height = 6, width = 8)

write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_240221.txt"), append = T)
with(import, coxph(OS ~ PC_M_60 + strata(PROTOCOL_REV) )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_240221.txt"), append = T, col_names = T)
with(import, coxph(PFS ~ PC_M_60 + strata(PROTOCOL_REV) )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_240221.txt"), append = T, col_names = T)

#______ LDH __________ 

gg <- ggsurvplot(survfit(OS ~ import$LDH_UL, data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("LDH_UL=0", "LDH_UL=1"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = print(gg), filename = "LDH_UL_OS.png", path = outpath, dpi = 300, height = 6, width = 8)
gg <- ggsurvplot(survfit(PFS ~ import$LDH_UL, data = import) , pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("LDH_UL=0", "LDH_UL=1"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = print(gg), filename = "LDH_UL_PFS.png", path = outpath, dpi = 300, height = 6, width = 8)

write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_240221.txt"), append = T)
with(import, coxph(OS ~ LDH_UL + strata(PROTOCOL_REV) )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_240221.txt"), append = T, col_names = T)
with(import, coxph(PFS ~ LDH_UL + strata(PROTOCOL_REV) )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_240221.txt"), append = T, col_names = T)

#______ ISS __________ 

gg <- ggsurvplot(survfit(OS ~ import$ISS, data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("ISS 1", "ISS 2", "ISS 3"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = print(gg), filename = "ISS_OS.png", path = outpath, dpi = 300, height = 6, width = 8)
gg <- ggsurvplot(survfit(PFS ~ import$ISS, data = import) , pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("ISS 1", "ISS 2", "ISS 3"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = print(gg), filename = "ISS_PFS.png", path = outpath, dpi = 300, height = 6, width = 8)

write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_240221.txt"), append = T)
with(import, coxph(OS ~ ISS + strata(PROTOCOL_REV) )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_240221.txt"), append = T, col_names = T)
with(import, coxph(PFS ~ ISS + strata(PROTOCOL_REV) )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_240221.txt"), append = T, col_names = T)

#______ R-ISS __________ 

gg <- ggsurvplot(survfit(OS ~ import$R_ISS, data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("R_ISS 1", "R_ISS 2", "R_ISS 3"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = print(gg), filename = "R_ISS_OS.png", path = outpath, dpi = 300, height = 6, width = 8)
gg <- ggsurvplot(survfit(PFS ~ import$R_ISS, data = import) , pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("R_ISS 1", "R_ISS 2", "R_ISS 3"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = print(gg), filename = "R_ISS_PFS.png", path = outpath, dpi = 300, height = 6, width = 8)

write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_240221.txt"), append = T)
with(import, coxph(OS ~ R_ISS + strata(PROTOCOL_REV) )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_240221.txt"), append = T, col_names = T)
with(import, coxph(PFS ~ R_ISS + strata(PROTOCOL_REV) )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_240221.txt"), append = T, col_names = T)

#______ Ig_Isotype __________ 

gg <- ggsurvplot(survfit(OS ~ import$IG_ISOTYPE_REV, data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("BJ", "IgA", "IgD", "IgG", "IgM"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = print(gg), filename = "IG_ISOTYPE_REV_OS.png", path = outpath, dpi = 300, height = 6, width = 8)
gg <- ggsurvplot(survfit(PFS ~ import$IG_ISOTYPE_REV, data = import) , pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("BJ", "IgA", "IgD", "IgG", "IgM"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = print(gg), filename = "IG_ISOTYPE_REV_PFS.png", path = outpath, dpi = 300, height = 6, width = 8)

write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_240221.txt"), append = T)
with(import, coxph(OS ~ IG_ISOTYPE_REV + strata(PROTOCOL_REV) )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_240221.txt"), append = T, col_names = T)
with(import, coxph(PFS ~ IG_ISOTYPE_REV + strata(PROTOCOL_REV) )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_240221.txt"), append = T, col_names = T)

#______ Light-Chain __________ 

gg <- ggsurvplot(survfit(OS ~ import$LIGHT_CHAIN, data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("Kappa Light Chain", "Lambda Light Chain"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = print(gg), filename = "LIGHT_CHAIN_OS.png", path = outpath, dpi = 300, height = 6, width = 8)
gg <- ggsurvplot(survfit(PFS ~ import$LIGHT_CHAIN, data = import) , pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("Kappa Light Chain", "Lambda Light Chain"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = print(gg), filename = "LIGHT_CHAIN_PFS.png", path = outpath, dpi = 300, height = 6, width = 8)


write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_240221.txt"), append = T)
with(import, coxph(OS ~ LIGHT_CHAIN + strata(PROTOCOL_REV) )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_240221.txt"), append = T, col_names = T)
with(import, coxph(PFS ~ LIGHT_CHAIN + strata(PROTOCOL_REV) )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_240221.txt"), append = T, col_names = T)


#______ TX yes/no ________ 

gg <- ggsurvplot(survfit(OS ~ import$TX_I_LINE_1_0, data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("TX_I_LINE= 0", "TX_I_LINE= 1"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = print(gg), filename = "TX_I_LINE_1_0_OS.png", path = outpath, dpi = 300, height = 6, width = 8)
gg <- ggsurvplot(survfit(PFS ~ import$TX_I_LINE_1_0, data = import) , pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("TX_I_LINE= 0", "TX_I_LINE= 1"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = print(gg), filename = "TX_I_LINE_1_0_PFS.png", path = outpath, dpi = 300, height = 6, width = 8)

write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_240221.txt"), append = T)
with(import, coxph(OS ~ TX_I_LINE_1_0 + strata(PROTOCOL_REV) )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_240221.txt"), append = T, col_names = T)
with(import, coxph(PFS ~ TX_I_LINE_1_0 + strata(PROTOCOL_REV) )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_240221.txt"), append = T, col_names = T)

#______ number TX __________ 

gg <- ggsurvplot(survfit(OS ~ import$NUMBERS_TX_FRONT_LINE, data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("NO Front Line ASCT", "DOUBLE Front Line ASCT", "SINGLE Front Line ASCT"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = print(gg), filename = "NUMBERS_TX_FRONT_LINE_OS.png", path = outpath, dpi = 300, height = 6, width = 8)
gg <- ggsurvplot(survfit(PFS ~ import$NUMBERS_TX_FRONT_LINE, data = import) , pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("NO Front Line ASCT", "DOUBLE Front Line ASCT", "SINGLE Front Line ASCT"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = print(gg), filename = "NUMBERS_TX_FRONT_LINE_PFS.png", path = outpath, dpi = 300, height = 6, width = 8)


write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_240221.txt"), append = T)
with(import, coxph(OS ~ NUMBERS_TX_FRONT_LINE + strata(PROTOCOL_REV) )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_240221.txt"), append = T, col_names = T)
with(import, coxph(PFS ~ NUMBERS_TX_FRONT_LINE + strata(PROTOCOL_REV) )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_240221.txt"), append = T, col_names = T)

#_______ t(4;14) _______ 

gg <- ggsurvplot(survfit(OS ~ import$FISH_T_4_14, data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("FISH_T_4_14=0", "FISH_T_4_14=1" ), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = print(gg), filename = "FISH_T_4_14_OS.png", path = outpath, dpi = 300, height = 6, width = 8)
gg <- ggsurvplot(survfit(PFS ~ import$FISH_T_4_14, data = import) , pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("FISH_T_4_14=0", "FISH_T_4_14=1"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = print(gg), filename = "FISH_T_4_14_PFS.png", path = outpath, dpi = 300, height = 6, width = 8)

write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_240221.txt"), append = T)
with(import, coxph(OS ~ FISH_T_4_14 + strata(PROTOCOL_REV) )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_240221.txt"), append = T, col_names = T)
with(import, coxph(PFS ~ FISH_T_4_14 + strata(PROTOCOL_REV) )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_240221.txt"), append = T, col_names = T)

#_______ t(11;14) _______ 

gg <- ggsurvplot(survfit(OS ~ import$FISH_T_11_14, data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("FISH_T_11_14=0", "FISH_T_11_14=1" ), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = print(gg), filename = "FISH_T_11_14_OS.png", path = outpath, dpi = 300, height = 6, width = 8)
gg <- ggsurvplot(survfit(PFS ~ import$FISH_T_11_14, data = import) , pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("FISH_T_11_14=0", "FISH_T_11_14=1"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = print(gg), filename = "FISH_T_11_14_PFS.png", path = outpath, dpi = 300, height = 6, width = 8)

write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_240221.txt"), append = T)
with(import, coxph(OS ~ FISH_T_11_14 + strata(PROTOCOL_REV) )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_240221.txt"), append = T, col_names = T)
with(import, coxph(PFS ~ FISH_T_11_14 + strata(PROTOCOL_REV) )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_240221.txt"), append = T, col_names = T)

#_______ t(14;16) _______ 

gg <- ggsurvplot(survfit(OS ~ import$FISH_T_14_16, data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("FISH_T_14_16=0", "FISH_T_14_16=1" ), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = print(gg), filename = "FISH_T_14_16_OS.png", path = outpath, dpi = 300, height = 6, width = 8)
gg <- ggsurvplot(survfit(PFS ~ import$FISH_T_14_16, data = import) , pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("FISH_T_14_16=0", "FISH_T_14_16=1"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = print(gg), filename = "FISH_T_14_16_PFS.png", path = outpath, dpi = 300, height = 6, width = 8)

write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_240221.txt"), append = T)
with(import, coxph(OS ~ FISH_T_14_16 + strata(PROTOCOL_REV) )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_240221.txt"), append = T, col_names = T)
with(import, coxph(PFS ~ FISH_T_14_16 + strata(PROTOCOL_REV) )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_240221.txt"), append = T, col_names = T)

#_______ t(14;20) _______ 

gg <- ggsurvplot(survfit(OS ~ import$FISH_T_14_20, data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("FISH_T_14_20=0", "FISH_T_14_20=1" ), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = print(gg), filename = "FISH_T_14_20_OS.png", path = outpath, dpi = 300, height = 6, width = 8)
gg <- ggsurvplot(survfit(PFS ~ import$FISH_T_14_20, data = import) , pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("FISH_T_14_20=0", "FISH_T_14_20=1"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = print(gg), filename = "FISH_T_14_20_PFS.png", path = outpath, dpi = 300, height = 6, width = 8)

write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_240221.txt"), append = T)
with(import, coxph(OS ~ FISH_T_14_20 + strata(PROTOCOL_REV) )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_240221.txt"), append = T, col_names = T)
with(import, coxph(PFS ~ FISH_T_14_20 + strata(PROTOCOL_REV) )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_240221.txt"), append = T, col_names = T)

#_______ t(6;14) _______ 

gg <- ggsurvplot(survfit(OS ~ import$FISH_T_6_14, data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("FISH_T_6_14=0", "FISH_T_6_14=1" ), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = print(gg), filename = "FISH_T_6_14_OS.png", path = outpath, dpi = 300, height = 6, width = 8)
gg <- ggsurvplot(survfit(PFS ~ import$FISH_T_6_14, data = import) , pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("FISH_T_6_14=0", "FISH_T_6_14=1"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = print(gg), filename = "FISH_T_6_14_PFS.png", path = outpath, dpi = 300, height = 6, width = 8)

write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_240221.txt"), append = T)
with(import, coxph(OS ~ FISH_T_6_14 + strata(PROTOCOL_REV) )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_240221.txt"), append = T, col_names = T)
with(import, coxph(PFS ~ FISH_T_6_14 + strata(PROTOCOL_REV) )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_240221.txt"), append = T, col_names = T)

#_______ Hyperdiploidy _______ 

gg <- ggsurvplot(survfit(OS ~ import$HyperDiploidy, data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("HyperDiploidy=0", "HyperDiploidy=1" ), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = print(gg), filename = "HyperDiploidy_OS.png", path = outpath, dpi = 300, height = 6, width = 8)
gg <- ggsurvplot(survfit(PFS ~ import$HyperDiploidy, data = import) , pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("HyperDiploidy=0", "HyperDiploidy=1"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = print(gg), filename = "HyperDiploidy_PFS.png", path = outpath, dpi = 300, height = 6, width = 8)


write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_240221.txt"), append = T)
with(import, coxph(OS ~ HyperDiploidy + strata(PROTOCOL_REV) )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_240221.txt"), append = T, col_names = T)
with(import, coxph(PFS ~ HyperDiploidy + strata(PROTOCOL_REV) )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_240221.txt"), append = T, col_names = T)


#______ DEL 17p _______ 

gg <- ggsurvplot(survfit(OS ~ import$"DEL_maj-broad_chr_17p", data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("Chr 17p Delection= 0", "Chr 17p Delection= 1" ), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = print(gg), filename = "DEL_maj-broad_chr_17p_OS.png", path = outpath, dpi = 300, height = 6, width = 8)
gg <- ggsurvplot(survfit(PFS ~ import$"DEL_maj-broad_chr_17p", data = import) , pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("Chr 17p Delection= 0", "Chr 17p Delection= 1"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = print(gg), filename = "DEL_maj-broad_chr_17p_PFS.png", path = outpath, dpi = 300, height = 6, width = 8)

write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_240221.txt"), append = T)
with(import, coxph(OS ~ `DEL_maj-broad_chr_17p` + strata(PROTOCOL_REV) )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_240221.txt"), append = T, col_names = T)
with(import, coxph(PFS ~ `DEL_maj-broad_chr_17p` + strata(PROTOCOL_REV) )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_240221.txt"), append = T, col_names = T)

#______ DEL TP53 ________ 

gg <- ggsurvplot(survfit(OS ~ import$"DEL_maj-focal_TP53", data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("TP53 Delection= 0", "TP53 Delection= 1" ), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = print(gg), filename = "DEL_maj-focal_TP53_OS.png", path = outpath, dpi = 300, height = 6, width = 8)
gg <- ggsurvplot(survfit(PFS ~ import$"DEL_maj-focal_TP53", data = import) , pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("TP53 Delection= 0", "TP53 Delection= 1"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = print(gg), filename = "DEL_maj-focal_TP53_PFS.png", path = outpath, dpi = 300, height = 6, width = 8)

write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_240221.txt"), append = T)
with(import, coxph(OS ~ `DEL_maj-focal_TP53` + strata(PROTOCOL_REV) )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_240221.txt"), append = T, col_names = T)
with(import, coxph(PFS ~ `DEL_maj-focal_TP53` + strata(PROTOCOL_REV) )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_240221.txt"), append = T, col_names = T)

#______ AMP 1q ________ 

gg <- ggsurvplot(survfit(OS ~ import$"AMP_maj-broad_chr_1q", data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("Chr 1q Broad_Amplification= 0", "Chr 1q Broad_Amplification= 1" ), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = print(gg), filename = "AMP_maj-broad_chr_1_OS.png", path = outpath, dpi = 300, height = 6, width = 8)
gg <- ggsurvplot(survfit(PFS ~ import$"AMP_maj-broad_chr_1q", data = import) , pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("Chr 1q Broad_Amplification= 0", "Chr 1q Broad_Amplification= 1"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = print(gg), filename = "AMP_maj-broad_chr_1_PFS.png", path = outpath, dpi = 300, height = 6, width = 8)

write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_240221.txt"), append = T)
with(import, coxph(OS ~ `AMP_maj-broad_chr_1q` + strata(PROTOCOL_REV) )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_240221.txt"), append = T, col_names = T)
with(import, coxph(PFS ~ `AMP_maj-broad_chr_1q` + strata(PROTOCOL_REV) )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_240221.txt"), append = T, col_names = T)



# 1q MMrisk focal+broad 

gg <- ggsurvplot(survfit(OS ~ import$MMrisk_1q_all, data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("Chr1q Broad&Focal Amp=0", "Chr1q Broad&Focal Amp=1" ), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = print(gg), filename = "MMrisk_1q_all_OS.png", path = outpath, dpi = 300, height = 6, width = 8)
gg <- ggsurvplot(survfit(PFS ~ import$MMrisk_1q_all, data = import) , pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("Chr1q Broad&Focal Amp=0", "Chr1q Broad&Focal Amp=1"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = print(gg), filename = "MMrisk_1q_all_PFS.png", path = outpath, dpi = 300, height = 6, width = 8)

write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_240221.txt"), append = T)
with(import, coxph(OS ~ MMrisk_1q_all + strata(PROTOCOL_REV) )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_240221.txt"), append = T, col_names = T)
with(import, coxph(PFS ~ MMrisk_1q_all + strata(PROTOCOL_REV) )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_240221.txt"), append = T, col_names = T)


#______ DEL 13 ________ 

gg <- ggsurvplot(survfit(OS ~ import$"DEL_maj-broad_chr_13q", data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("Chr 13q Broad Delection= 0", "Chr 13q Broad Delection= 1" ), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = print(gg), filename = "DEL_maj-broad_chr_13q_OS.png", path = outpath, dpi = 300, height = 6, width = 8)
gg <- ggsurvplot(survfit(PFS ~ import$"DEL_maj-broad_chr_13q", data = import) , pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("Chr 13q Broad Delection= 0", "Chr 13q Broad Delection= 1"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = print(gg), filename = "DEL_maj-broad_chr_13q_PFS.png", path = outpath, dpi = 300, height = 6, width = 8)


write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_240221.txt"), append = T)
with(import, coxph(OS ~ `DEL_maj-broad_chr_13q` + strata(PROTOCOL_REV) )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_240221.txt"), append = T, col_names = T)
with(import, coxph(PFS ~ `DEL_maj-broad_chr_13q` + strata(PROTOCOL_REV) )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_240221.txt"), append = T, col_names = T)


# 13 MMrisk focal+broad

gg <- ggsurvplot(survfit(OS ~ import$MMrisk_13_all, data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("Chr 13q Broad&Focal Del= 0", "Chr 13q Broad&Focal Del= 1" ), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = print(gg), filename = "MMrisk_13_all_OS.png", path = outpath, dpi = 300, height = 6, width = 8)
gg <- ggsurvplot(survfit(PFS ~ import$MMrisk_13_all, data = import) , pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("Chr 13q Broad&Focal Del= 0", "Chr 13q Broad&Focal Del= 1"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = print(gg), filename = "MMrisk_13_all_PFS.png", path = outpath, dpi = 300, height = 6, width = 8)

write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_240221.txt"), append = T)
with(import, coxph(OS ~ MMrisk_13_all + strata(PROTOCOL_REV) )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_240221.txt"), append = T, col_names = T)
with(import, coxph(PFS ~ MMrisk_13_all + strata(PROTOCOL_REV) )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_240221.txt"), append = T, col_names = T)

#______ MM RISK ________ 

gg <- ggsurvplot(survfit(OS ~ import$MMrisk_CLASSIFIER_ALL, data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("1q&13+", "1q/13+", "1q&13-" ), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = print(gg), filename = "MMrisk_CLASSIFIER_ALL_OS.png", path = outpath, dpi = 300, height = 6, width = 8)
gg <- ggsurvplot(survfit(PFS ~ import$MMrisk_CLASSIFIER_ALL, data = import) , pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("1q&13+", "1q/13+","1q/13-" ), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = print(gg), filename = "MMrisk_CLASSIFIER_ALL_PFS.png", path = outpath, dpi = 300, height = 6, width = 8)


write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_240221.txt"), append = T)
with(import, coxph(OS ~ MMrisk_CLASSIFIER_ALL + strata(PROTOCOL_REV) )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_240221.txt"), append = T, col_names = T)
with(import, coxph(PFS ~ MMrisk_CLASSIFIER_ALL + strata(PROTOCOL_REV) )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_240221.txt"), append = T, col_names = T)



#______ Risk 1 ______ 

gg <- ggsurvplot(survfit(OS ~ import$MMrisk1, data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("1q&13+ =0", "1q&13+ =1" ), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = print(gg), filename = "MMrisk1_OS.png", path = outpath, dpi = 300, height = 6, width = 8)
gg <- ggsurvplot(survfit(PFS ~ import$MMrisk1, data = import) , pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("1q&13+ =0", "1q&13+ =1"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = print(gg), filename = "MMrisk1_PFS.png", path = outpath, dpi = 300, height = 6, width = 8)

write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_240221.txt"), append = T)
with(import, coxph(OS ~ MMrisk1 + strata(PROTOCOL_REV) )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_240221.txt"), append = T, col_names = T)
with(import, coxph(PFS ~ MMrisk1 + strata(PROTOCOL_REV) )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_240221.txt"), append = T, col_names = T)

#_______ Risk 2 ______ 

gg <- ggsurvplot(survfit(OS ~ import$MMrisk2, data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("1q/13+ =0", "1q/13+ =1" ), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = print(gg), filename = "MMrisk2_OS.png", path = outpath, dpi = 300, height = 6, width = 8)
gg <- ggsurvplot(survfit(PFS ~ import$MMrisk2, data = import) , pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("1q/13+ =0", "1q/13+ =1"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = print(gg), filename = "MMrisk2_PFS.png", path = outpath, dpi = 300, height = 6, width = 8)


write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_240221.txt"), append = T)
with(import, coxph(OS ~ MMrisk2 + strata(PROTOCOL_REV) )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_240221.txt"), append = T, col_names = T)
with(import, coxph(PFS ~ MMrisk2 + strata(PROTOCOL_REV) )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_240221.txt"), append = T, col_names = T)


#_______ Risk 3 _______ 

gg <- ggsurvplot(survfit(OS ~ import$MMrisk3, data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("1q&13- =0", "1q&13- =1" ), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = print(gg), filename = "MMrisk3_OS.png", path = outpath, dpi = 300, height = 6, width = 8)
gg <- ggsurvplot(survfit(PFS ~ import$MMrisk3, data = import) , pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("1q&13- =0", "1q&13- =1"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = print(gg), filename = "MMrisk3_PFS.png", path = outpath, dpi = 300, height = 6, width = 8)


write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_240221.txt"), append = T)
with(import, coxph(OS ~ MMrisk3 + strata(PROTOCOL_REV) )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_240221.txt"), append = T, col_names = T)
with(import, coxph(PFS ~ MMrisk3 + strata(PROTOCOL_REV) )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_240221.txt"), append = T, col_names = T)


#____ MM RISK a/b _______

gg <- ggsurvplot(survfit(OS ~ import$MMrisk_AB_ALL, data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("1q&13+", "only 13 Del", "only 1q Amp", "1q&13-" ), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = print(gg), filename = "MMrisk_AB_ALL_OS.png", path = outpath, dpi = 300, height = 6, width = 8)
gg <- ggsurvplot(survfit(PFS ~ import$MMrisk_AB_ALL, data = import) , pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("1q&13+", "only 13 Del", "only 1q Amp", "1q&13-"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = print(gg), filename = "MMrisk_AB_ALL_PFS.png", path = outpath, dpi = 300, height = 6, width = 8)



write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_240221.txt"), append = T)
with(import, coxph(OS ~ MMrisk_AB_ALL + strata(PROTOCOL_REV) )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_240221.txt"), append = T, col_names = T)
with(import, coxph(PFS ~ MMrisk_AB_ALL + strata(PROTOCOL_REV) )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_240221.txt"), append = T, col_names = T)



#____ Risk_1 vs only_1q _______

gg <- ggsurvplot(survfit(OS ~ import$risk1_VS_1q_VS_other, data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("1q&13+", "only 1q Amp", "other" ), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = print(gg), filename = "risk1_VS_1q_VS_other_OS.png", path = outpath, dpi = 300, height = 6, width = 8)
gg <- ggsurvplot(survfit(PFS ~ import$risk1_VS_1q_VS_other, data = import) , pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("1q&13+", "only 1q amp", "other"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = print(gg), filename = "risk1_VS_1q_VS_other_PFS.png", path = outpath, dpi = 300, height = 6, width = 8)



write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_240221.txt"), append = T)
with(import, coxph(OS ~ risk1_VS_1q_VS_other + strata(PROTOCOL_REV) )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_240221.txt"), append = T, col_names = T)
with(import, coxph(PFS ~ risk1_VS_1q_VS_other + strata(PROTOCOL_REV) )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_240221.txt"), append = T, col_names = T)



imp1q <- filter(import, import$"AMP_maj-broad_chr_1q" == 1)

OS1q <- Surv(imp1q$OS_MONTHS, imp1q$OS_EVENT)
PFS1q <- Surv(imp1q$PFS_I_MONTHS, imp1q$PFS_I_EVENT)

gg <- ggsurvplot(survfit(OS1q ~ imp1q$risk1_VS_1q_VS_other, data = imp1q) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("1q&13+", "only 1q Amp" ), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = print(gg), filename = "risk1_VS_1q_VS_other_OS1q.png", path = outpath, dpi = 300, height = 6, width = 8)
gg <- ggsurvplot(survfit(PFS1q ~ imp1q$risk1_VS_1q_VS_other, data = imp1q) , pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("1q&13+", "only 1q Amp"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = print(gg), filename = "risk1_VS_1q_VS_other_PFS1q.png", path = outpath, dpi = 300, height = 6, width = 8)


write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_240221.txt"), append = T)
with(import, coxph(OS1q ~ risk1_VS_1q_VS_other + strata(PROTOCOL_REV) )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_240221.txt"), append = T, col_names = T)
with(import, coxph(PFS1q ~ risk1_VS_1q_VS_other + strata(PROTOCOL_REV) )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_240221.txt"), append = T, col_names = T)


#______ ultra MMrisk classificator ______

gg <- ggsurvplot(survfit(OS ~ import$ULTRA_MMrisk, data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("ULTRA_MMrisk=1", "ULTRA_MMrisk=2", "ULTRA_MMrisk=3" ), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = print(gg), filename = "ULTRA_MMrisk_OS.png", path = outpath, dpi = 300, height = 6, width = 8)
gg <- ggsurvplot(survfit(PFS ~ import$ULTRA_MMrisk, data = import), pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("ULTRA_MMrisk=1", "ULTRA_MMrisk=2", "ULTRA_MMrisk=3"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = print(gg), filename = "ULTRA_MMrisk_PFS.png", path = outpath, dpi = 300, height = 6, width = 8)

write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_240221.txt"), append = T)
with(import, coxph(OS ~ ULTRA_MMrisk + strata(PROTOCOL_REV) )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_240221.txt"), append = T, col_names = T)
with(import, coxph(PFS ~ ULTRA_MMrisk + strata(PROTOCOL_REV) )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_240221.txt"), append = T, col_names = T)



#______ ultra high  ______

gg <- ggsurvplot(survfit(OS ~ import$ULTRA_High, data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("ULTRA_High=0", "ULTRA_High=1"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = print(gg), filename = "ULTRA_High_OS.png", path = outpath, dpi = 300, height = 6, width = 8)
gg <- ggsurvplot(survfit(PFS ~ import$ULTRA_High, data = import), pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("ULTRA_High=0", "ULTRA_High=1"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = print(gg), filename = "ULTRA_High_PFS.png", path = outpath, dpi = 300, height = 6, width = 8)

write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_240221.txt"), append = T)
with(import, coxph(OS ~ ULTRA_High + strata(PROTOCOL_REV) )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_240221.txt"), append = T, col_names = T)
with(import, coxph(PFS ~ ULTRA_High + strata(PROTOCOL_REV) )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_240221.txt"), append = T, col_names = T)



#______ ultra low  ______

gg <- ggsurvplot(survfit(OS ~ import$ULTRA_Low, data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("ULTRA_Low=0", "ULTRA_Low=1"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = print(gg), filename = "ULTRA_High_OS.png", path = outpath, dpi = 300, height = 6, width = 8)
gg <- ggsurvplot(survfit(PFS ~ import$ULTRA_Low, data = import), pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("ULTRA_Low=0", "ULTRA_Low=1"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)
ggsave(plot = print(gg), filename = "ULTRA_High_PFS.png", path = outpath, dpi = 300, height = 6, width = 8)

write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_240221.txt"), append = T)
with(import, coxph(OS ~ ULTRA_Low + strata(PROTOCOL_REV) )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_240221.txt"), append = T, col_names = T)
with(import, coxph(PFS ~ ULTRA_Low + strata(PROTOCOL_REV) )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_240221.txt"), append = T, col_names = T)




#=======================================================================
#========================== MULTIVARIATE ===============================
#=======================================================================


#_____ multivariate ______
write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_240221.txt"), append = T)
write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_240221.txt"), append = T)
write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_240221.txt"), append = T)
write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_240221.txt"), append = T)

write_tsv(data.frame("MULTIVARIATA 1 "),paste0(outpath,"report_cox_analysis_240221.txt"), append = T)
with(import, coxph(OS ~ OLD_0_1 + B2M + ALBUMIN + HB_m_105 + PLT_m_150 + MMrisk1 + TX_I_LINE_1_0 + `DEL_maj-focal_TP53` + FISH_T_4_14 + FISH_T_14_16 + FISH_T_14_20 + `DEL_maj-broad_chr_1p` + HyperDiploidy + strata(PROTOCOL_REV) )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_240221.txt"), append = T, col_names = T)
with(import, coxph(PFS ~OLD_0_1 + B2M + ALBUMIN + HB_m_105 + PLT_m_150 + MMrisk1 + TX_I_LINE_1_0 + `DEL_maj-focal_TP53` + FISH_T_4_14 + FISH_T_14_16 + FISH_T_14_20 + `DEL_maj-broad_chr_1p` + HyperDiploidy + strata(PROTOCOL_REV) )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_240221.txt"), append = T, col_names = T)


write_tsv(data.frame("MULTIVARIATA 2 "),paste0(outpath,"report_cox_analysis_240221.txt"), append = T)
with(import, coxph(OS ~ B2M + PLT_m_150 + MMrisk1 + TX_I_LINE_1_0 + `DEL_maj-focal_TP53` + FISH_T_4_14 +  FISH_T_14_20 + `DEL_maj-broad_chr_1p` + strata(PROTOCOL_REV) )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_240221.txt"), append = T, col_names = T)
with(import, coxph(PFS ~ B2M + PLT_m_150 + MMrisk1 + TX_I_LINE_1_0 +  FISH_T_14_20 + strata(PROTOCOL_REV) )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_240221.txt"), append = T, col_names = T)

write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_240221.txt"), append = T)

write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_240221.txt"), append = T)
write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_240221.txt"), append = T)


write_tsv(data.frame("MULTIVARIATA ULTRA HIGH 1 "),paste0(outpath,"report_cox_analysis_240221.txt"), append = T)
with(import, coxph(OS ~ OLD_0_1 +  HB_m_105 + PLT_m_150 + ULTRA_High + TX_I_LINE_1_0 + `DEL_maj-focal_TP53` + FISH_T_4_14 + FISH_T_14_16 + FISH_T_14_20 + `DEL_maj-broad_chr_1p` + HyperDiploidy + strata(PROTOCOL_REV) )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_240221.txt"), append = T, col_names = T)
with(import, coxph(PFS ~OLD_0_1 +  HB_m_105 + PLT_m_150 + ULTRA_High + TX_I_LINE_1_0 + `DEL_maj-focal_TP53` + FISH_T_4_14 + FISH_T_14_16 + FISH_T_14_20 + `DEL_maj-broad_chr_1p` + HyperDiploidy + strata(PROTOCOL_REV) )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_240221.txt"), append = T, col_names = T)

write_tsv(data.frame("\n"),paste0(outpath,"report_cox_analysis_240221.txt"), append = T)

write_tsv(data.frame("MULTIVARIATA ULTRA HIGH 2 "),paste0(outpath,"report_cox_analysis_240221.txt"), append = T)
with(import, coxph(OS ~ PLT_m_150 + ULTRA_High + TX_I_LINE_1_0 + `DEL_maj-focal_TP53` + FISH_T_4_14 + FISH_T_14_20 + `DEL_maj-broad_chr_1p` + strata(PROTOCOL_REV) )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_240221.txt"), append = T, col_names = T)
with(import, coxph(PFS ~ PLT_m_150 + ULTRA_High + TX_I_LINE_1_0 + FISH_T_14_20 + strata(PROTOCOL_REV) )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)%>% write_tsv(paste0(outpath,"report_cox_analysis_240221.txt"), append = T, col_names = T)



