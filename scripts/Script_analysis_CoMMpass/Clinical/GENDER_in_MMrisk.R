library(data.table)
library(tidyverse)
library(survival)
library(survminer)

import <- fread("D:/analisi_in_corso/risk_1q_13/complete_database_1q_13_181019.txt") 
import$PROTOCOL %>% table

#______ MM RISK CREATION _____

import$AMP_1q_genes_all <- with(import, 
                                ifelse(`AMP_maj_focal_ANP32E`==1 | `AMP_maj_focal_MCL1` ==1 | `AMP_maj_focal_CKS1B`==1, 1,0 ))


import$MMrisk_1q_all <- ifelse( import$`AMP_maj_broad_chr_1q` == 1 | import$AMP_1q_genes_all ==1, 1, 0)
import$MMrisk_13_all <- ifelse( import$`DEL_maj_broad_chr_13q` ==1 | import$`DEL_maj_focal_RB1` ==1, 1, 0)

import$MMrisk_CLASSIFIER_ALL <- 3  - import$MMrisk_1q_all - import$MMrisk_13_all

#_______ PFS OS creation ______

OS <- Surv(import$OS_months, import$OS_event_death)
PFS <- Surv(import$PFS_I_months, import$PFS_I_event)


################# GENDER ANLYSIS ####################

ggsurvplot(survfit(OS ~ Sex, data = import), pval = T, risk.table = T)
ggsurvplot(survfit(PFS ~ Sex, data = import), pval = T, risk.table = T)


#____ sex in risk groups ____
import$R3_M_F <- ifelse(import$MMrisk_CLASSIFIER_ALL==3 & import$Sex == "F", "risk3_F",
                        ifelse(import$MMrisk_CLASSIFIER_ALL==3 & import$Sex == "M", "risk3_M", paste0("risk_",import$MMrisk_CLASSIFIER_ALL)))

import$R3_M_F %>% table

ggsurvplot(survfit(OS ~ R3_M_F, data = import), pval = T, risk.table = T)
ggsurvplot(survfit(PFS ~ R3_M_F, data = import), pval = T, risk.table = T)


coxph(OS ~ import$DEL_maj_broad_chr_17p,data = import) %>% summary
coxph(PFS ~ import$DEL_maj_broad_chr_17p,data = import) %>% summary


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

#_________ MM risk per GENDER_______

fun <- Vectorize( function(x,y) {if(x=="M") {
  paste0("risk_",y,"_M")
} else {
  paste0("risk_",y,"_F")
}
})

import$MMrisk_GENDER <- fun(import$Sex, import$MMrisk_CLASSIFIER_ALL)

import$MMrisk_GENDER %>% table

ggsurvplot(survfit(OS ~ MMrisk_GENDER, data = import), pval = T, risk.table = T)
ggsurvplot(survfit(PFS ~ MMrisk_GENDER, data = import), pval = T, risk.table = T)

#_________ MM risk per GENDER per ISS_______

fun2 <- Vectorize( function(x,y,z) {if(x=="M") {
  paste0("risk_",y,"_M_",z)
} else {
  paste0("risk_",y,"_F_",z)
}
})

import$MMrisk_GENDER_ISS <- fun2(import$Sex, import$MMrisk_CLASSIFIER_ALL, import$ISS)
import$MMrisk_GENDER_ISS %>% table
import$MMrisk_GENDER_ISS %>% unique

import$MMrisk_GENDER_ISS[grep("_NA",import$MMrisk_GENDER_ISS)] <- NA

ggsurvplot(survfit(OS ~ MMrisk_GENDER_ISS, data = import), pval = T, risk.table = T)
ggsurvplot(survfit(PFS ~ MMrisk_GENDER_ISS, data = import), pval = T, risk.table = T)


#================== Classificator mmrisk / GENDER / ISS =======================



import$CLASS <- with(import, ifelse(MMrisk_GENDER_ISS == "risk_3_F_1", "lowest",
                                  ifelse(MMrisk_GENDER_ISS == "risk_1_M_3" | MMrisk_GENDER_ISS == "risk_1_F_3", "highest", "mid"))
                     )
import$CLASS = relevel(import$CLASS %>% as.factor, ref = "mid")


ggsurvplot(survfit(OS ~ CLASS, data = import), pval = T, risk.table = T)
ggsurvplot(survfit(PFS ~ CLASS, data = import), pval = T, risk.table = T)

coxph(OS ~ import$CLASS,data = import) %>% summary
coxph(PFS ~ import$CLASS,data = import) %>% summary



lowest <- import %>% filter(CLASS=="lowest")
highest <- import %>% filter(CLASS=="highest")


lowest$DEL_maj_focal_TP53 %>% table
highest$DEL_maj_focal_TP53 %>% table

table(import$DEL_maj_focal_TP53, import$CLASS)

gmodels::CrossTable(import$DEL_maj_focal_TP53, import$CLASS, fisher = T)


import$PLT %>% median(na.rm = T)
import$PLT_high %>% median(na.rm = T)

###################################################################################
#############################    COMMPASS    ######################################
###################################################################################


import <- fread("D:/analisi_in_corso/risk_1q_13/CoMMpass/complete_database_1q_13_CoMMpass.txt") 


#______ MM RISK CREATION _____

import$AMP_1q_genes_all <- with(import, 
                                ifelse(`AMP_maj_focal_ANP32E`==1 | `AMP_maj_focal_MCL1` ==1 | `AMP_maj_focal_CKS1B`==1, 1,0 ))


import$MMrisk_1q_all <- ifelse( import$`AMP_maj_broad_chr_1q` == 1 | import$AMP_1q_genes_all ==1, 1, 0)
import$MMrisk_13_all <- ifelse( import$`DEL_maj_broad_chr_13q` ==1 | import$`DEL_maj_focal_RB1` ==1, 1, 0)

import$MMrisk_CLASSIFIER_ALL <- 3  - import$MMrisk_1q_all - import$MMrisk_13_all

# add sex and clincal data
import_clinical <- fread("D:/CoMMpass_data/IA13a/Clinical_Flat_Files/CoMMpass_IA13_FlatFiles/MMRF_CoMMpass_IA13_PER_PATIENT.csv")

import_clinical$Study_Visit_iD <- paste0(import_clinical$PUBLIC_ID, "_1_BM")

import2 <- left_join(import, import_clinical,  by = "Study_Visit_iD")

import <- import2

table(import$MMrisk_CLASSIFIER_ALL, import$D_PT_gender)

import$Sex <- ifelse(import$D_PT_gender==1, "M", "F")


OS <- Surv((import$ttcos/30.5) %>% ceiling, import$censos)
PFS <- Surv((import$ttcpfs/30.5) %>% ceiling, import$censpfs)


ggsurvplot(survfit(OS ~ Sex, data = import), pval = T, risk.table = T)
ggsurvplot(survfit(PFS ~ Sex, data = import), pval = T, risk.table = T)


#_________ MM risk per GENDER_______

fun <- Vectorize( function(x,y) {if(x=="M") {
  paste0("risk_",y,"_M")
} else {
  paste0("risk_",y,"_F")
}
})

import$MMrisk_GENDER <- fun(import$Sex, import$MMrisk_CLASSIFIER_ALL)

import$MMrisk_GENDER %>% table

ggsurvplot(survfit(OS ~ MMrisk_GENDER, data = import), pval = T, risk.table = T)
ggsurvplot(survfit(PFS ~ MMrisk_GENDER, data = import), pval = T, risk.table = T)

#_________ MM risk per GENDER per ISS_______

fun2 <- Vectorize( function(x,y,z) {if(x=="M") {
  paste0("risk_",y,"_M_",z)
} else {
  paste0("risk_",y,"_F_",z)
}
})

import$MMrisk_GENDER_ISS <- fun2(import$Sex, import$MMrisk_CLASSIFIER_ALL, import$D_PT_iss)

import$MMrisk_GENDER_ISS %>% table
import$MMrisk_GENDER_ISS %>% unique

import$MMrisk_GENDER_ISS[grep("_NA",import$MMrisk_GENDER_ISS)] <- NA

ggsurvplot(survfit(OS ~ MMrisk_GENDER_ISS, data = import), pval = T, risk.table = T)
ggsurvplot(survfit(PFS ~ MMrisk_GENDER_ISS, data = import), pval = T, risk.table = T)


#================== Classificator mmrisk / GENDER / ISS =======================



import$CLASS <- with(import, ifelse(MMrisk_GENDER_ISS == "risk_3_F_1", "lowest",
                                    ifelse(MMrisk_GENDER_ISS == "risk_1_M_3" | MMrisk_GENDER_ISS == "risk_1_F_3", "highest", "mid"))
)
import$CLASS = relevel(import$CLASS %>% as.factor, ref = "mid")


ggsurvplot(survfit(OS ~ CLASS, data = import), pval = T, risk.table = T)
ggsurvplot(survfit(PFS ~ CLASS, data = import), pval = T, risk.table = T)

coxph(OS ~ import$CLASS,data = import) %>% summary
coxph(PFS ~ import$CLASS,data = import) %>% summary


