library(data.table)
library(tidyverse)


# import gistic arm calls (in integer CN)

GISTIC_armCN <- fread("data/GISTIC_run132274/validrun1.broad_values_by_arm.txt")

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


# assign gain, high-copy gain and no gain labels to each event based on thresholds

# gain > 0.5 amplified copies
# high copy gain > 1.5 amplified copies 
# no gain <= 0.5 amplified copies

plotdata_melt$status <- ifelse(plotdata_melt$value>0.5 & plotdata_melt$value<= 1.5, "gain",
                               ifelse(plotdata_melt$value>1.5, "high_copy_gain", "no_gain"))

chr1q <- plotdata_melt %>% filter(variable=="chr_1q")

table(chr1q$status)


clinical_bolo <- fread("results/complete_database_BO_1q_13.txt")


merge <- left_join(chr1q, clinical_bolo, by=c("sample"="SNP_FILE_NAME"))

#============================== surv analysis on 1q gain VS high copy gain =========================
library(survival)
library(survminer)

import <- merge
import$PFS_I_MONTHS

OS <- Surv(import$OS_MONTHS, import$OS_EVENT)
PFS <- Surv(import$PFS_I_MONTHS, import$PFS_I_EVENT)

# example: TP53 deletion
ggsurvplot(survfit(OS ~ import$DEL_maj_focal_TP53, data = import), pval = T, risk.table = T)
ggsurvplot(survfit(PFS ~ import$DEL_maj_focal_TP53, data = import), pval = T, risk.table = T)

coxph(OS ~ import$DEL_maj_focal_TP53,data = import) %>% summary
coxph(PFS ~ import$DEL_maj_focal_TP53,data = import) %>% summary


#=========== 1q gain VS high copy gain analysis =============== 
ggsurvplot(survfit(OS ~ import$status, data = import), pval = T, risk.table = T) + xlab("OS")
ggsurvplot(survfit(PFS ~ import$status, data = import), pval = T, risk.table = T) + xlab("PFS")

coxph(OS ~ import$status,data = import) %>% summary
coxph(PFS ~ import$status,data = import) %>% summary


#___________________ 1q only______________
importAMP1q <- import %>% filter(status != "no_gain")

OS_1q <- Surv(importAMP1q$OS_MONTHS, importAMP1q$OS_EVENT)
PFS_1q <- Surv(importAMP1q$PFS_I_MONTHS, importAMP1q$PFS_I_EVENT)

ggsurvplot(survfit(OS_1q ~ importAMP1q$status, data = importAMP1q), pval = T, risk.table = T) + xlab("OS")
ggsurvplot(survfit(PFS_1q ~ importAMP1q$status, data = importAMP1q), pval = T, risk.table = T) + xlab("PFS")

coxph(OS_1q ~ importAMP1q$status,data = importAMP1q) %>% summary
coxph(PFS_1q ~ importAMP1q$status,data = importAMP1q) %>% summary

# not significant difference both in PFS and OS


