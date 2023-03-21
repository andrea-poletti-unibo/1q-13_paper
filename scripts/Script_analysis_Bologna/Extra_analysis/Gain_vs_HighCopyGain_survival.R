library(data.table)
library(tidyverse)


# import gistic arm calls (in integer CN)

# GISTIC_armLogR <- fread("D:/analisi_in_corso/risk_1q_13/GISTIC/127469_run2_RawCentered/run2.broad_values_by_arm.txt")
# GISTIC_armCN <- fread("D:/analisi_in_corso/risk_1q_13/old_batches/GISTIC/127469_run2_RawCentered/run2.broad_values_by_arm.txt")

GISTIC_armCN <- fread("D:/analisi_in_corso/risk_1q_13/GISTIC/132274_validrun1/validrun1.broad_values_by_arm.txt")


# transpose and rename column
armsCN <- t(as.matrix(GISTIC_armCN, rownames=1)) %>% as.data.frame
names(armsCN) <- paste0("chr_", names(armsCN))


#======= explorative analysis on clonality =========

plotdata <- armsCN
plotdata$sample <- rownames(plotdata) %>% str_replace("_\\.CytoScanHD_Array\\.","") %>% str_replace("_\\.GenomeWideSNP_6\\.","")
plotdata_melt <- melt(plotdata, id.vars= "sample")


############# function to find all the peaks ##############
dens <- density(plotdata_melt$value)
plot(dens)
maximums<-c()
for(i in 2:(length(dens$y)-1)){
  if(dens$y[i]>dens$y[i-1] & dens$y[i]>dens$y[i+1]){
    maximums<-c(maximums,i)
  }
}
abline(v=dens$x[maximums],lty=2)
dens$x[maximums] # 1.108 for del and 1.136 for amps
###########################################################

# ~ 1.12 is the scaling factor: mean between amp and del peak

########### SCALE the signal for the derived scaling factor #############
armsCN_scaled <- apply(armsCN, 2, function(x) x/1.12 ) %>% as.data.frame()

# re-melt the data matrix
plotdata <- armsCN_scaled

plotdata$sample <- rownames(plotdata) %>% str_replace("^X","") %>% str_replace("\\.CytoScanHD_Array\\.","(CytoScanHD_Array)") %>% str_replace("\\.GenomeWideSNP_6\\.","(GenomeWideSNP_6)")  %>% str_replace("\\.1035_"," 1035_")   

plotdata_melt <- melt(plotdata, id.vars= "sample")


# assign clonal, subclonal major and minor labels to each event based on threshold
plotdata_melt$status <- ifelse(plotdata_melt$value>0.5 & plotdata_melt$value<1.5, "gain",
                               ifelse(plotdata_melt$value>1.5, "high_copy_gain", "no_gain"))

chr1q <- plotdata_melt %>% filter(variable=="chr_1q")

table(chr1q$status)


clinical_bolo <- fread("D:/analisi_in_corso/risk_1q_13/complete_database_1q_13_181019.txt")


merge <- left_join(chr1q, clinical_bolo, by=c("sample"="nome_CEL_sample"))

#============================== surv analysis on 1q gain VS high copy gain =========================
library(survival)
library(survminer)

import <- merge

OS <- Surv(import$OS_months, import$OS_event_death)
PFS <- Surv(import$PFS_I_months, import$PFS_I_event)

ggsurvplot(survfit(OS ~ import$DEL_maj_focal_TP53, data = import), pval = T, risk.table = T)
ggsurvplot(survfit(PFS ~ import$DEL_maj_focal_TP53, data = import), pval = T, risk.table = T)

coxph(OS ~ import$DEL_maj_focal_TP53,data = import) %>% summary
coxph(PFS ~ import$DEL_maj_focal_TP53,data = import) %>% summary



ggsurvplot(survfit(OS ~ import$status, data = import), pval = T, risk.table = T) + xlab("OS")
ggsurvplot(survfit(PFS ~ import$status, data = import), pval = T, risk.table = T) + xlab("PFS")

coxph(OS ~ import$status,data = import) %>% summary
coxph(PFS ~ import$status,data = import) %>% summary


#___________________ 1q only______________
importAMP1q <- import %>% filter(status != "no_gain")

OS_1q <- Surv(importAMP1q$OS_months, importAMP1q$OS_event_death)
PFS_1q <- Surv(importAMP1q$PFS_I_months, importAMP1q$PFS_I_event)

ggsurvplot(survfit(OS_1q ~ importAMP1q$status, data = importAMP1q), pval = T, risk.table = T) + xlab("OS")
ggsurvplot(survfit(PFS_1q ~ importAMP1q$status, data = importAMP1q), pval = T, risk.table = T) + xlab("PFS")

coxph(OS_1q ~ importAMP1q$status,data = importAMP1q) %>% summary
coxph(PFS_1q ~ importAMP1q$status,data = importAMP1q) %>% summary
