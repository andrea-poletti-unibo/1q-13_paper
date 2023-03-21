
library(data.table)
library(tidyverse)


# import gistic arm calls (in integer CN)
# GISTIC_armLogR <- fread("D:/analisi_in_corso/risk_1q_13/GISTIC/127469_run2_RawCentered/run2.broad_values_by_arm.txt")
# GISTIC_armCN <- fread("D:/analisi_in_corso/risk_1q_13/old_batches/GISTIC/127469_run2_RawCentered/run2.broad_values_by_arm.txt")

GISTIC_armCN <- fread("D:/analisi_in_corso/risk_1q_13/CoMMpass/175814_50_validrunCoMMpass/jointsize50.broad_values_by_arm.txt")

# transpose and rename column
armsCN <- t(as.matrix(GISTIC_armCN, rownames=1)) %>% as.data.frame
names(armsCN) <- paste0("chr_", names(armsCN))


#======= explorative analysis on clonality =========
plotdata <- armsCN
plotdata$sample <- rownames(plotdata) 
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

dens$x[maximums] # peaks at -0.954 for clonal dels and 0.912 for clonal amps
###########################################################

(0.954 + 0.912)/ 2
# 0.933 is the scaling factor: mean between amp and del peak

########### SCALE the signal for the scaling factor #############
armsCN_scaled <- apply(armsCN, 2, function(x) x/0.933 ) %>% as.data.frame()


# re-melt the data matrix
plotdata <- armsCN_scaled
plotdata$sample <- rownames(plotdata) %>% str_replace("_\\.CytoScanHD_Array\\.","") %>% str_replace("_\\.GenomeWideSNP_6\\.","")
plotdata_melt <- melt(plotdata, id.vars= "sample")

# assign clonal, subclonal major and minor labels to each event based on threshold
plotdata_melt$status <- ifelse(abs(plotdata_melt$value)>=0.9, 
                               "clonal",
                               ifelse(abs(plotdata_melt$value)<= 0.1,
                                      "normal",
                                      ifelse(abs(plotdata_melt$value)>= 0.5,
                                             "subclonal_major","subclonal_minor")))


#=========== TABLES and PLOTS ===========
#_____ table ______
table_data <- plotdata_melt

table_data <- filter(table_data, status != "normal")

table_data$event <- ifelse(table_data$value > 0.1, "AMP", 
                              ifelse(table_data$value < -0.1, "DEL", "no_CNA"))

summary_table <- table_data %>% group_by(variable, event) %>% 
  summarise(clonal= sum(status=="clonal"), 
            subclonal_major = sum(status=="subclonal_major"), 
            subclonal_minor = sum(status=="subclonal_minor"), 
            clonal_perc= sum(status=="clonal")/514, 
            subclonal_major_perc = sum(status=="subclonal_major")/514,
            subclonal_minor_perc = sum(status=="subclonal_minor")/514
            ) %>% 
  mutate(MAJOR_TOT= clonal+subclonal_major,
         MAJOR_TOT_PERC=MAJOR_TOT/514 )

# write_tsv(summary_table, "D:/analisi_in_corso/risk_1q_13/GISTIC/summarytable_clonal-subclonal_calls.txt")

plotdata_melt %>% ggplot(aes(value)) + geom_density() # scaled around CN1 for del peak and CN3 for amp peak

#============== plots for each event in every sample by arm ==============

17400 # end of chr 10
# plotdata_melt[17400,]
# plotdata_melt[17401,]

# use 17400 to cut the dataframe in two parts, and plot them separately

# plot 2.a
plotdata_melt[1:17400,] %>% filter(status != "normal") %>% ggplot(aes(variable, value)) + 
  geom_violin(scale="count")+
  geom_jitter(aes(colour= status), width=0.2, alpha=0.5, size=1.2) +
  geom_hline(yintercept = 0.1, linetype=1) +
  geom_hline(yintercept = -0.1, linetype=1) + ylim(-2,6) + geom_hline(yintercept = 0.9, linetype=3) + 
  geom_hline(yintercept = -0.9, linetype=3) + geom_hline(yintercept = -0.5, linetype=3) + geom_hline(yintercept = 0.5, linetype=3)

# plot 2.b
plotdata_melt[17400:33930,] %>% filter(status != "normal") %>% ggplot(aes(variable, value)) + 
  geom_violin(scale="count") +
  geom_jitter(aes(colour= status), width=0.2, alpha=0.5, size=1.2) +
  geom_hline(yintercept = 0.1, linetype=1) +
  geom_hline(yintercept = -0.1, linetype=1) + ylim(-2,6) + geom_hline(yintercept = 0.9, linetype=3) + 
  geom_hline(yintercept = -0.9, linetype=3) + geom_hline(yintercept = -0.5, linetype=3) + geom_hline(yintercept = 0.5, linetype=3)

# all chrs
plotdata_melt %>% filter(status != "normal") %>% ggplot(aes(variable, value + 2)) + 
  geom_violin(scale="count")+
  labs(x="chromosome arm", y=" Copy Number") +
  geom_jitter(aes(colour= status), width=0.2, alpha=0.5, size=1.2) +
  geom_hline(yintercept = 2.1, linetype=1) +
  geom_hline(yintercept = 1.9, linetype=1) + ylim(0.5,5) + 
  geom_hline(yintercept = 2.9, linetype=3) + 
  geom_hline(yintercept = 1.1, linetype=3) + 
  geom_hline(yintercept = 1.5, linetype=3) + 
  geom_hline(yintercept = 2.5, linetype=3) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust= 1))

ggsave(filename = "violin_plot_all_broad_calls_CoMMpass.pdf",path = "C:/Users/mm_gr/Alma Mater Studiorum Università di Bologna/PROJECT 1q & 13 - Documenti/PAPER_FIGURES_OFFICIAL/", 
       device = "pdf",units = "in",dpi = 320, width = 16, height = 8)





###################### BROAD CALLS ##########################


#============== CALLS 90% (CLONAL) ===================
# call a clonal alteration if > 0.9 from diploid 

CALLS_AMP_CLONAL <- apply(as.matrix(armsCN_scaled), c(1,2), function(x) if(x>0.9) 1 else 0) %>% as.data.frame()
names(CALLS_AMP_CLONAL) <- paste0("AMP_clon_broad_", names(CALLS_AMP_CLONAL))

library(gplots)
heatmap.2(as.matrix(CALLS_AMP_CLONAL), col="redgreen", trace= "none")

write.table(CALLS_AMP_CLONAL,"D:/analisi_in_corso/risk_1q_13/GISTIC/BROAD_CALLS/AMP_CALLS_clonal_broad.txt", sep = "\t", quote = F, row.names = T)


CALLS_DEL_CLONAL <- apply(as.matrix(armsCN_scaled), c(1,2), function(x) if(x< -0.9) 1 else 0) %>% as.data.frame()
names(CALLS_DEL_CLONAL) <- paste0("DEL_clon_broad_", names(CALLS_DEL_CLONAL))

library(gplots)
heatmap.2(as.matrix(CALLS_DEL_CLONAL), col="redgreen", trace= "none")

write.table(CALLS_DEL_CLONAL,"D:/analisi_in_corso/risk_1q_13/GISTIC/BROAD_CALLS/DEL_CALLS_clonal_broad.txt", sep = "\t", quote = F, row.names = T)

#============== CALLS 50% (MAJOR CLONE) ===================
# call a clonal alteration if > 0.5 from diploid 

CALLS_AMP_MAJOR <- apply(as.matrix(armsCN_scaled), c(1,2), function(x) if(x>0.5) 1 else 0) %>% as.data.frame()
names(CALLS_AMP_MAJOR) <- paste0("AMP_maj_broad_", names(CALLS_AMP_MAJOR))

library(gplots)
heatmap.2(as.matrix(CALLS_AMP_MAJOR), col="redgreen", trace= "none")

write.table(CALLS_AMP_MAJOR,"D:/analisi_in_corso/risk_1q_13/GISTIC/BROAD_CALLS/AMP_CALLS_major_broad.txt", sep = "\t", quote = F, row.names = T)


CALLS_DEL_MAJOR <- apply(as.matrix(armsCN_scaled), c(1,2), function(x) if(x< -0.5) 1 else 0) %>% as.data.frame()
names(CALLS_DEL_MAJOR) <- paste0("DEL_maj_broad_", names(CALLS_DEL_MAJOR))

library(gplots)
heatmap.2(as.matrix(CALLS_DEL_MAJOR), col="redgreen", trace= "none")

write.table(CALLS_DEL_MAJOR,"D:/analisi_in_corso/risk_1q_13/GISTIC/BROAD_CALLS/DEL_CALLS_major_broad.txt", sep = "\t", quote = F, row.names = T)


#============== CALLS 10% (ALL) ===================
# call a clonal alteration if > 0.1 from diploid 

CALLS_AMP_ALL <- apply(as.matrix(armsCN_scaled), c(1,2), function(x) if(x>0.1) 1 else 0) %>% as.data.frame()
names(CALLS_AMP_ALL) <- paste0("AMP_all_broad_", names(CALLS_AMP_ALL))

library(gplots)
heatmap.2(as.matrix(CALLS_AMP_ALL), col="redgreen", trace= "none")

write.table(CALLS_AMP_ALL,"D:/analisi_in_corso/risk_1q_13/GISTIC/BROAD_CALLS/AMP_CALLS_all_broad.txt", sep = "\t", quote = F, row.names = T)


CALLS_DEL_ALL <- apply(as.matrix(armsCN_scaled), c(1,2), function(x) if(x< -0.1) 1 else 0) %>% as.data.frame()
names(CALLS_DEL_ALL) <- paste0("DEL_all_broad_", names(CALLS_DEL_ALL))

library(gplots)
heatmap.2(as.matrix(CALLS_DEL_ALL), col="redgreen", trace= "none")

write.table(CALLS_DEL_ALL,"D:/analisi_in_corso/risk_1q_13/GISTIC/BROAD_CALLS/DEL_CALLS_all_broad.txt", sep = "\t", quote = F, row.names = T)



# # Cluster complete
# MAT <- cbind(CALLS_AMP_ALL, CALLS_DEL_ALL)
# 
# heatmap.2(as.matrix(MAT), col=c("black","green"), trace= "none")
# 
# MAT$SUM <- apply(as.matrix(MAT),1,sum)
# 
# MAT_0 <- MAT %>% filter(SUM==0)
