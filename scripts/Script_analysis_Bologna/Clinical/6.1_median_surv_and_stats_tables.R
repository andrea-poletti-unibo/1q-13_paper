library(data.table)
library(tidyverse)
library(survival)
library(survminer)
library(RODBC)
library(broom)

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

import$MMrisk_AB_ALL <- ifelse(import$MMrisk_CLASS == 2 & import$MMrisk_1q_all ==1,
                               "2_1q", 
                               ifelse(import$MMrisk_CLASS == 2 & import$MMrisk_13_all ==1,
                                      "2_13",
                                      import$MMrisk_CLASS))

import$risk1_VS_1q_VS_other <- recode(import$MMrisk_AB_ALL, "1"="1q&13", "2_1q"="1q_only", "2_13"="other", "3"="other")


import$MMrisk1 <- ifelse(import$MMrisk_CLASS==1, 1, 0)
import$MMrisk2 <- ifelse(import$MMrisk_CLASS==2, 1, 0)
import$MMrisk3 <- ifelse(import$MMrisk_CLASS==3, 1, 0)



import$ULTRA_MMrisk <- ifelse(import$ISS == 3 & import$MMrisk_CLASS == 1,
                              "Ultra_High",
                              ifelse(import$ISS==1 & import$MMrisk_CLASS == 3,
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

import$Induction_therapy <- ifelse(!import$INDUCTION_THERAPY_FRONT_LINE %in% c("VCD", "VTD", "TD" ), "other", import$INDUCTION_THERAPY_FRONT_LINE %>% as.character() )
import$Induction_therapy %>% table

import$Induction_response_M_VGPR <- ifelse(import$INDUCTION_RESPONSE %in% c("VGPR", "nCR", "CR", "sCR"), 1,0)
import$Induction_response_M_CR <- ifelse(import$INDUCTION_RESPONSE %in% c("nCR", "CR", "sCR"), 1,0)

################################################################################ 
################################ SURV ANALYSIS #################################
################################################################################ 


outpath <- "C:/Users/andre/Alma Mater Studiorum Università di Bologna/PROJECT 1q & 13 - Documenti/Survival_analysis/Bologna_dataset/analysis_240321/"
dir.create(outpath)

write_tsv(data.frame("median_months_curves"),paste0(outpath,"median_months_curves_290321.txt"), append = T)



#_______ PROTOCOL ________ 

survfit(OS ~ import$PROTOCOL_REV, data = import) %>% summary %>% .$table %>% as.data.frame() %>% cbind(surv="OS", .) %>% rownames_to_column(var = "var") %>% write_tsv(paste0(outpath,"median_months_curves_290321.txt"), append = T, col_names = T)
gg <- ggsurvplot(survfit(OS ~ import$PROTOCOL_REV, data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("AMB", "BO2005", "EMN02"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)

survfit(PFS ~ import$PROTOCOL_REV, data = import) %>% summary %>% .$table %>% as.data.frame() %>% cbind(surv="PFS", .) %>% rownames_to_column(var = "var") %>% write_tsv(paste0(outpath,"median_months_curves_290321.txt"), append = T, col_names = T)
gg <- ggsurvplot(survfit(PFS ~ import$PROTOCOL_REV, data = import) , pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("AMB", "BO2005", "EMN02"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)


#______ MM RISK ________ 

survfit(OS ~ import$MMrisk_class, data = import) %>% summary %>% .$table %>% as.data.frame() %>% cbind(surv="OS", .) %>% rownames_to_column(var = "var") %>% write_tsv(paste0(outpath,"median_months_curves_290321.txt"), append = T, col_names = T)
gg <- ggsurvplot(survfit(OS ~ import$MMrisk_class, data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("1q&13+", "1q/13+", "1q&13-" ), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)


survfit(PFS ~ import$MMrisk_class, data = import) %>% summary %>% .$table %>% as.data.frame() %>% cbind(surv="PFS", .) %>% rownames_to_column(var = "var") %>% write_tsv(paste0(outpath,"median_months_curves_290321.txt"), append = T, col_names = T)
gg <- ggsurvplot(survfit(PFS ~ import$MMrisk_class, data = import) , pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("1q&13+", "1q/13+","1q&13-" ), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)



#______ ultra MMrisk classificator ______

survfit(OS ~ import$ULTRA_MMrisk, data = import) %>% summary %>% .$table %>% as.data.frame() %>% cbind(surv="OS", .) %>% rownames_to_column(var = "var") %>% write_tsv(paste0(outpath,"median_months_curves_290321.txt"), append = T, col_names = T)
gg <- ggsurvplot(survfit(OS ~ import$ULTRA_MMrisk, data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("ULTRA_MMrisk=1", "ULTRA_MMrisk=2", "ULTRA_MMrisk=3" ), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)

survfit(PFS ~ import$ULTRA_MMrisk, data = import) %>% summary %>% .$table %>% as.data.frame() %>% cbind(surv="PFS", .) %>% rownames_to_column(var = "var") %>% write_tsv(paste0(outpath,"median_months_curves_290321.txt"), append = T, col_names = T)
gg <- ggsurvplot(survfit(PFS ~ import$ULTRA_MMrisk, data = import), pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("ULTRA_MMrisk=1", "ULTRA_MMrisk=2", "ULTRA_MMrisk=3"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)




#______ Risk 1 ______ 

survfit(OS ~ import$MMrisk1, data = import) %>% summary %>% .$table %>% as.data.frame() %>% cbind(surv="OS", .) %>% rownames_to_column(var = "var") %>% write_tsv(paste0(outpath,"median_months_curves_290321.txt"), append = T, col_names = T)
gg <- ggsurvplot(survfit(OS ~ import$MMrisk1, data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("1q&13+ =0", "1q&13+ =1" ), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)

survfit(PFS ~ import$MMrisk1, data = import) %>% summary %>% .$table %>% as.data.frame() %>% cbind(surv="PFS", .) %>% rownames_to_column(var = "var") %>% write_tsv(paste0(outpath,"median_months_curves_290321.txt"), append = T, col_names = T)
gg <- ggsurvplot(survfit(PFS ~ import$MMrisk1, data = import) , pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("1q&13+ =0", "1q&13+ =1"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)


#_______ Risk 2 ______ 

survfit(OS ~ import$MMrisk2, data = import) %>% summary %>% .$table %>% as.data.frame() %>% cbind(surv="OS", .) %>% rownames_to_column(var = "var") %>% write_tsv(paste0(outpath,"median_months_curves_290321.txt"), append = T, col_names = T)
gg <- ggsurvplot(survfit(OS ~ import$MMrisk2, data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("1q/13+ =0", "1q/13+ =1" ), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)

survfit(PFS ~ import$MMrisk2, data = import) %>% summary %>% .$table %>% as.data.frame() %>% cbind(surv="PFS", .) %>% rownames_to_column(var = "var") %>% write_tsv(paste0(outpath,"median_months_curves_290321.txt"), append = T, col_names = T)
gg <- ggsurvplot(survfit(PFS ~ import$MMrisk2, data = import) , pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("1q/13+ =0", "1q/13+ =1"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)


#_______ Risk 3 _______ 

survfit(OS ~ import$MMrisk3, data = import) %>% summary %>% .$table %>% as.data.frame() %>% cbind(surv="OS", .) %>% rownames_to_column(var = "var") %>% write_tsv(paste0(outpath,"median_months_curves_290321.txt"), append = T, col_names = T)
gg <- ggsurvplot(survfit(OS ~ import$MMrisk3, data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("1q&13- =0", "1q&13- =1" ), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)

survfit(PFS ~ import$MMrisk3, data = import) %>% summary %>% .$table %>% as.data.frame() %>% cbind(surv="PFS", .) %>% rownames_to_column(var = "var") %>% write_tsv(paste0(outpath,"median_months_curves_290321.txt"), append = T, col_names = T)
gg <- ggsurvplot(survfit(PFS ~ import$MMrisk3, data = import) , pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("1q&13- =0", "1q&13- =1"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)

