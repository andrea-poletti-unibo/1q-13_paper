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

outpath <- "C:/Users/andre/Alma Mater Studiorum Università di Bologna/PROJECT 1q & 13 - Documenti/Survival_analysis/CoMMpass_dataset/analysis_230321/"
dir.create(outpath)

write_tsv(data.frame("median_months_curves"),paste0(outpath,"median_months_curves_commpass_070421.txt"), append = T)

#______ MM RISK ________ 

survfit(OS ~ import$MMrisk_class, data = import) %>% summary %>% .$table %>% as.data.frame() %>% cbind(surv="OS", .) %>% rownames_to_column(var = "var") %>% write_tsv(paste0(outpath,"median_months_curves_commpass_070421.txt"), append = T, col_names = T)
gg <- ggsurvplot(survfit(OS ~ import$MMrisk_class, data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("1q&13+", "1q/13+", "1q&13-" ), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)


survfit(PFS ~ import$MMrisk_class, data = import) %>% summary %>% .$table %>% as.data.frame() %>% cbind(surv="PFS", .) %>% rownames_to_column(var = "var") %>% write_tsv(paste0(outpath,"median_months_curves_commpass_070421.txt"), append = T, col_names = T)
gg <- ggsurvplot(survfit(PFS ~ import$MMrisk_class, data = import) , pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("1q&13+", "1q/13+","1q&13-" ), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)



#______ ultra MMrisk classificator ______

survfit(OS ~ import$ULTRA_MMrisk, data = import) %>% summary %>% .$table %>% as.data.frame() %>% cbind(surv="OS", .) %>% rownames_to_column(var = "var") %>% write_tsv(paste0(outpath,"median_months_curves_commpass_070421.txt"), append = T, col_names = T)
gg <- ggsurvplot(survfit(OS ~ import$ULTRA_MMrisk, data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("ULTRA_MMrisk=1", "ULTRA_MMrisk=2", "ULTRA_MMrisk=3" ), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)

survfit(PFS ~ import$ULTRA_MMrisk, data = import) %>% summary %>% .$table %>% as.data.frame() %>% cbind(surv="PFS", .) %>% rownames_to_column(var = "var") %>% write_tsv(paste0(outpath,"median_months_curves_commpass_070421.txt"), append = T, col_names = T)
gg <- ggsurvplot(survfit(PFS ~ import$ULTRA_MMrisk, data = import), pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("ULTRA_MMrisk=1", "ULTRA_MMrisk=2", "ULTRA_MMrisk=3"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)




#______ Risk 1 ______ 

survfit(OS ~ import$MMrisk1, data = import) %>% summary %>% .$table %>% as.data.frame() %>% cbind(surv="OS", .) %>% rownames_to_column(var = "var") %>% write_tsv(paste0(outpath,"median_months_curves_commpass_070421.txt"), append = T, col_names = T)
gg <- ggsurvplot(survfit(OS ~ import$MMrisk1, data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("1q&13+ =0", "1q&13+ =1" ), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)

survfit(PFS ~ import$MMrisk1, data = import) %>% summary %>% .$table %>% as.data.frame() %>% cbind(surv="PFS", .) %>% rownames_to_column(var = "var") %>% write_tsv(paste0(outpath,"median_months_curves_commpass_070421.txt"), append = T, col_names = T)
gg <- ggsurvplot(survfit(PFS ~ import$MMrisk1, data = import) , pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("1q&13+ =0", "1q&13+ =1"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)


#_______ Risk 2 ______ 

survfit(OS ~ import$MMrisk2, data = import) %>% summary %>% .$table %>% as.data.frame() %>% cbind(surv="OS", .) %>% rownames_to_column(var = "var") %>% write_tsv(paste0(outpath,"median_months_curves_commpass_070421.txt"), append = T, col_names = T)
gg <- ggsurvplot(survfit(OS ~ import$MMrisk2, data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("1q/13+ =0", "1q/13+ =1" ), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)

survfit(PFS ~ import$MMrisk2, data = import) %>% summary %>% .$table %>% as.data.frame() %>% cbind(surv="PFS", .) %>% rownames_to_column(var = "var") %>% write_tsv(paste0(outpath,"median_months_curves_commpass_070421.txt"), append = T, col_names = T)
gg <- ggsurvplot(survfit(PFS ~ import$MMrisk2, data = import) , pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("1q/13+ =0", "1q/13+ =1"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)


#_______ Risk 3 _______ 

survfit(OS ~ import$MMrisk3, data = import) %>% summary %>% .$table %>% as.data.frame() %>% cbind(surv="OS", .) %>% rownames_to_column(var = "var") %>% write_tsv(paste0(outpath,"median_months_curves_commpass_070421.txt"), append = T, col_names = T)
gg <- ggsurvplot(survfit(OS ~ import$MMrisk3, data = import) , pval = T, risk.table = T, xlab = "OS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("1q&13- =0", "1q&13- =1" ), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)

survfit(PFS ~ import$MMrisk3, data = import) %>% summary %>% .$table %>% as.data.frame() %>% cbind(surv="PFS", .) %>% rownames_to_column(var = "var") %>% write_tsv(paste0(outpath,"median_months_curves_commpass_070421.txt"), append = T, col_names = T)
gg <- ggsurvplot(survfit(PFS ~ import$MMrisk3, data = import) , pval = T, risk.table = T, xlab = "PFS", surv.median.line = "hv", break.time.by = 12, legend= c(0.9,0.9), legend.title="", legend.labs=c("1q&13- =0", "1q&13- =1"), tables.y.text = F, risk.table.y.text.col = TRUE, font.legend=c("bold"))
print(gg)

