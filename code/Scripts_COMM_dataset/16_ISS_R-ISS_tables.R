library(tidyverse)
library(data.table)
library(tableone)
library(kableExtra)


#______import clinical data already created
import <- data.table::fread("results/complete_database_CoMMpass_1q_13.txt") %>% as.data.frame()

outpath <- "results/Clinical_anlysis/CoMM_dataset/Multivariate/"
dir.create(outpath, recursive = T, showWarnings = F)

#_____ADDITIONAL DATA MANAGEMENT____
import$LDH_level<-ifelse(import$LDH_level=="",NA, import$LDH_level)
import$LDH_level<-as.factor(import$LDH_level)
levels(import$LDH_level)

import$R_ISS<-as.factor(import$R_ISS)
levels(import$R_ISS)

import$SEX <- ifelse(import$D_PT_gender==1, "M", "F") %>% as.factor()

import$SEX<-as.factor(import$SEX)
import$SEX<-relevel(import$SEX, ref = 2)

import$LDH <- import$D_LAB_chem_ldh

import$ASCT <- import$sctflag

import$FIRST_THERAPHY

import$ALBUMIN <- import$D_LAB_chem_albumin
import$ALBUMIN_m_3.5 <- ifelse(import$ALBUMIN < 35, 1,0)

import$B2M <- import$D_LAB_serum_beta2_microglobulin
import$B2M_level <- ifelse(import$B2M < 3.5, "<3.5", ifelse(import$B2M > 5.5, ">5.5", "normal")) %>% as.factor()
import$B2M_level <- import$B2M_level %>% relevel(ref = 3)
levels(import$B2M_level)

import$SeqFISH_Del_TP53 <- import$DEL_maj_focal_TP53

import$SeqFISH_Del_1p_CDKN2C <- ifelse(import$DEL_maj_focal_CDKN2C ==1, 1,0)


#============ create trasnlocations vars ============
import$IgH_translocation_type <-ifelse(import$SeqWGS_WHSC1_CALL==1, "t(4;14)",
                                       ifelse(import$SeqWGS_CCND1_CALL==1, "t(11;14)",
                                              ifelse(import$SeqWGS_CCND3_CALL==1,"t(6;14)",
                                                     ifelse(import$SeqWGS_MAF_CALL==1,"t(14;16)",
                                                            ifelse(import$SeqWGS_MAFB_CALL==1,"t(14;20)", 
                                                                   "no translocation")))))

table(import$IgH_translocation_type)
table(import$IgH_translocation_type)/840

import$IgH_translocation_type<-as.factor(import$IgH_translocation_type)
levels(import$IgH_translocation_type)

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

#_____ COMPUTE ISS ______

import$ISS <- ifelse(import$D_LAB_serum_beta2_microglobulin>= 5.5, 3,
                     ifelse(import$D_LAB_serum_beta2_microglobulin< 3.5 & import$D_LAB_chem_albumin >= 35, 1,2)) %>% as.factor() %>% relevel(ref = 2)
import$ISS
table(import$ISS)


import$R_ISS <- import$R_ISS %>% relevel(ref = 2)


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



###################################### TABLES ##########################################

data <- import

data$MMrisk <- data$MMrisk_CLASSIFIER_ALL %>% factor(levels = c("risk1","risk2","risk3"), labels = c("1q&13+","1q/13","1q&13-"))

data$R_ISS <- data$R_ISS %>% as.character()

data$MM_group %>% table

data$R_ISS %>% table(useNA = "always")
data$R_ISS[is.na(data$R_ISS)] <- "NA"

data$R_ISS_class <- data$R_ISS %>% factor(levels = c(1,2,3,"NA"), labels = c("R-ISS 1","R-ISS 2","R-ISS 3", "non-classified"))




data %>% select(MMrisk_class, MMrisk_class_t_CCND2, R_ISS_class) %>% 
    CreateTableOne(
        strata = "R_ISS_class", 
        vars = c("MMrisk_class", "MMrisk_class_t_CCND2"), 
        data = ., 
        test = TRUE, 
        addOverall = T, 
        includeNA = T) %>% 
    kableone() %>% 
    kable_classic() %>% kable_styling(full_width = F) %>% 
    save_kable(file = "results/Clinical_anlysis/CoMM_dataset/tables/R-ISS_and_1q13_table.html" )


data$ISS <- data$ISS %>% as.character()
data$ISS %>% table(useNA = "always")
data$ISS[is.na(data$ISS)] <- "NA"


data$ISS_class <- data$ISS %>% factor(levels = c(1,2,3,"NA"), labels = c("ISS 1","ISS 2","ISS 3", "non-classified"))


data %>% select(MMrisk_class, MMrisk_class_t_CCND2, ISS_class) %>% 
    CreateTableOne(
        strata = "ISS_class", 
        vars = c("MMrisk_class", "MMrisk_class_t_CCND2"), 
        data = ., 
        test = TRUE, 
        addOverall = T, 
        margin=1,
        includeNA = T) %>% 
    kableone() %>% 
    kable_classic() %>% kable_styling(full_width = F) %>% 
    save_kable(file = "results/Clinical_anlysis/CoMM_dataset/tables/ISS_and_1q13_table.html" )


tableone()


install.packages("Gmisc")
library(Gmisc)

mtcars %>% getDescriptionStatsBy(mpg, wt, am, by=am)

tab1 <- data %>% getDescriptionStatsBy(MMrisk_class , MMrisk_class_t_CCND2, by=ISS_class, header_count = TRUE, add_total_col = T, total_col_show_perc = F, hrzl_prop = T, percentage_sign = F, statistics = T, NEJMstyle=T) %>% htmlTable()  %>% addHtmlTableStyle(css.cell = c("width: 100;","width: 100;","width: 100;","width: 100;", "width: 100;","width: 100;")) 

tab1

sink("results/Clinical_anlysis/CoMM_dataset/tables/test.html")
print(tab1,type="html",useViewer=F)
sink()

