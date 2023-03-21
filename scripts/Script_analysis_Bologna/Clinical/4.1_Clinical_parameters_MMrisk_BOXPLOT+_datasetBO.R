
library(tidyverse)
library(data.table)
library(gridExtra)

data <- fread("C:/Users/emat/Alma Mater Studiorum Universit? di Bologna/PROJECT 1q & 13 - Documenti/complete_database_1q_13_181019.txt")

data$AMP_m

#____ Build CALLS for arms with broad+focal_____

#amps
data$AMP_maj_call_chr_1q <- ifelse( data$AMP_maj_broad_chr_1q == 1 | data$AMP_maj_focal_ANP32E== 1 | data$AMP_maj_focal_CKS1B ==1 | data$AMP_maj_focal_MCL1 ==1, 1, 0)

data$AMP_maj_call_chr_8q <- ifelse( data$AMP_maj_broad_chr_8q == 1 | data$AMP_maj_focal_MYC == 1, 1, 0)

#dels
data$DEL_maj_call_chr_17p <- ifelse( data$DEL_maj_broad_chr_17p == 1 | data$DEL_maj_focal_TP53 == 1, 1, 0)
data$DEL_maj_call_chr_1p  <- ifelse( data$DEL_maj_broad_chr_1p == 1 | data$DEL_maj_focal_CDKN2C == 1 | data$DEL_maj_focal_FAM46C == 1 | data$DEL_maj_focal_FAF1 == 1,  1, 0)
data$DEL_maj_call_chr_13q <- ifelse( data$DEL_maj_broad_chr_13q == 1 | data$DEL_maj_focal_RB1 == 1, 1, 0)
data$DEL_maj_call_chr_14q <- ifelse( data$DEL_maj_broad_chr_14q == 1 | data$DEL_maj_focal_TRAF3 == 1, 1, 0)
data$DEL_maj_call_chr_16q <- ifelse( data$DEL_maj_broad_chr_16q == 1 | data$DEL_maj_focal_CYLD == 1, 1, 0)

#_____ add group risk 1q & 13 _______
MM_risk <- ifelse( data$AMP_maj_call_chr_1q == 1 & data$DEL_maj_call_chr_13q ==1, 1,
                   ifelse(data$AMP_maj_call_chr_1q == 1 | data$DEL_maj_call_chr_13q ==1, 2, 3) ) 

data$MM_group_num <- MM_risk
data$MM_group <- MM_risk %>% factor(levels = c("1","2","3"), labels = c("1q&13+","1q/13","1q&13-"))

data$MM_risk_1 <- ifelse(MM_risk==1, 1,0)
data$MM_risk_2 <- ifelse(MM_risk==2, 1,0)
data$MM_risk_3 <- ifelse(MM_risk==3, 1,0)

#______ CLINICAL PARAMETERS BOXPLOTS _________

library(ggpubr)

my_comparisons <- list( c("1q&13-", "1q/13"), 
                        c("1q/13", "1q&13+"), 
                        c("1q&13-", "1q&13+") )

#______ albumin _______ ***

data$Albumin %>% summary

p1 <- data %>% ggplot(aes(MM_group, Albumin)) + 
  geom_violin() + 
  geom_boxplot(aes(fill=MM_group), alpha=0.5) + 
  stat_summary(geom="text", fun.y =quantile,
               aes(label=sprintf("%1.1f", ..y..), color=MM_group),
               position=position_nudge(y= 0.48), size=3) +
  stat_compare_means(comparisons = my_comparisons)+ 
  stat_compare_means(label.y.npc ="bottom") + xlab("") 


#______ calcium _______ **
data$Calcium %>% summary

p2 <- data %>% ggplot(aes(MM_group, Calcium)) + 
  geom_violin() + 
  geom_boxplot(aes(fill=MM_group), alpha=0.5) + 
  stat_summary(geom="text", fun.y = quantile,
               aes(label=sprintf("%1.1f", ..y..), color=MM_group),
               position=position_nudge(x=0.48), size=3) +
  stat_compare_means(comparisons = my_comparisons)+ 
  stat_compare_means(label.y.npc = "bottom") + xlab("") 


#______ creatinine _______ ns
data$Creatinine %>% summary

p3 <- data %>% ggplot(aes(MM_group, Creatinine)) + 
  geom_violin() + 
  geom_boxplot(aes(fill=MM_group), alpha=0.5) + 
  stat_summary(geom="text", fun.y = quantile,
               aes(label=sprintf("%1.1f", ..y..), color=MM_group),
               position=position_nudge(x=0.48), size=3) +
  stat_compare_means(comparisons = my_comparisons)+ 
  stat_compare_means(label.y.npc = "bottom") + xlab("") 


#______ hemoglobin _______ ***
data$HB %>% summary

p4 <- data %>% ggplot(aes(MM_group, HB)) + 
  geom_violin() + 
  geom_boxplot(aes(fill=MM_group), alpha=0.5) + 
  stat_summary(geom="text", fun.y = quantile,
               aes(label=sprintf("%1.1f", ..y..), color=MM_group),
               position=position_nudge(x=0.48), size=3) +
  stat_compare_means(comparisons = my_comparisons)+ 
  stat_compare_means(label.y.npc = "bottom") + xlab("") 


#______ LDH _______ ns
data$LDH %>% summary

p5 <- data %>% ggplot(aes(MM_group, LDH)) + 
  geom_violin() + 
  geom_boxplot(aes(fill=MM_group), alpha=0.5) + 
  stat_summary(geom="text", fun.y = quantile,
               aes(label=sprintf("%1.1f", ..y..), color=MM_group),
               position=position_nudge(x=0.48), size=3) +
  stat_compare_means(comparisons = my_comparisons)+ 
  stat_compare_means(label.y.npc = "bottom") + xlab("") 



#______ platlet _______ ***
data$PLTs %>% summary

p6 <- data %>% ggplot(aes(MM_group, PLTs)) + 
  geom_violin() + 
  geom_boxplot(aes(fill=MM_group), alpha=0.5) + 
  stat_summary(geom="text", fun.y = quantile,
               aes(label=sprintf("%1.1f", ..y..), color=MM_group),
               position=position_nudge(x=0.48), size=3) +
  stat_compare_means(comparisons = my_comparisons)+ 
  stat_compare_means(label.y.npc = "bottom") + xlab("") 




#______ beta2_microglobulin _______*** .
data$beta2m %>% summary

p7 <- data %>% ggplot(aes(MM_group, beta2m)) + 
  geom_violin() + 
  geom_boxplot(aes(fill=MM_group), alpha=0.5) + 
  stat_summary(geom="text", fun.y =quantile,
               aes(label=sprintf("%1.1f", ..y..), color=MM_group),
               position=position_nudge(x=0.48), size=3) +
  stat_compare_means(comparisons = my_comparisons)+ 
  stat_compare_means(label.y.npc = "bottom") + xlab("") 




Panel1 <- ggarrange(p1,p2,p3,p4,p5,p6,p7,
                    labels = c(paste0("a",1:7)),
                    ncol = 3, nrow = 3,
                    hjust = -0.5,
                    vjust = 1.5,
                    
)

Panel1 

ggsave("panel1.png", Panel1, path = "C:/Users/emat/Desktop/prova_clinica_bo", dpi = 300, width = 15, height = 21)


