library(data.table)
library(tidyverse)

import <- fread("focal_data_VS_june20.txt")

df <- import %>% select(`PAZ_001_E-197`:`PAZ_124_R-1267`) %>% t

colnames(df) <- import$relevant_genes

df <- df %>% as.data.frame()

df2 <- df %>% select(chr1q= `chr 1q22 (PDZK1, CD160, RNF115)`, chr13=`RB1 (LPAR6)`) %>% rownames_to_column(var = "sample")
df2$pt <- df2$sample %>% str_extract("PAZ_[0-9]+")

df_diagnos <- df2 %>% filter(grepl("_E", df2$sample))
df_relapse <- df2 %>% filter(grepl("_R", df2$sample))

df_all <- left_join(df_diagnos, df_relapse, by="pt",suffix = c("_diagnosis", "_relapse") )

def <- df_all %>% select(pt, chr1q_diagnosis,chr1q_relapse,chr13_diagnosis,chr13_relapse)

def$call_1q_diagnosis <- ifelse(def$chr1q_diagnosis> 2.5, "amp", ifelse(def$chr1q_diagnosis< 1.5, "del","normal"))
def$call_1q_relapse <- ifelse(def$chr1q_relapse> 2.5, "amp", ifelse(def$chr1q_relapse< 1.5, "del","normal"))

def$call_13_diagnosis <- ifelse(def$chr13_diagnosis> 2.5, "amp", ifelse(def$chr13_diagnosis< 1.5, "del","normal"))
def$call_13_relapse <- ifelse(def$chr13_relapse> 2.5, "amp", ifelse(def$chr13_relapse< 1.5, "del","normal"))

def$group_diagnosis <- ifelse(def$call_1q_diagnosis =="amp" & def$call_13_diagnosis =="del", "1q&13+", 
                              ifelse(def$call_1q_diagnosis =="normal" & def$call_13_diagnosis =="normal","1q&13-", "1q/13"))

def$group_relapse <- ifelse(def$call_1q_relapse =="amp" & def$call_13_relapse =="del", "1q&13+", 
                              ifelse(def$call_1q_relapse =="normal" & def$call_13_relapse =="normal","1q&13-", "1q/13"))



table(def$group_diagnosis, def$group_relapse)
