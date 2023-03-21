
library(tidyverse)
library(data.table)
library(gridExtra)

data <- fread("C:/Users/mm_gr//Alma Mater Studiorum Università di Bologna/PROJECT 1q & 13 - Documenti/complete_database_1q_13_CoMMpass_240121.txt")

#____ Build CALLS for arms with broad+focal_____

#amps
data$AMP_maj_call_chr_1q <- ifelse( data$AMP_maj_broad_chr_1q == 1 | data$AMP_maj_focal_ANP32E== 1 | data$AMP_maj_focal_CKS1B ==1 | data$AMP_maj_focal_MCL1 ==1, 1, 0)
data$AMP_maj_call_chr_8q <- ifelse( data$AMP_maj_broad_chr_8q == 1 | data$AMP_maj_focal_MYC == 1, 1, 0)

#dels
data$DEL_maj_call_chr_17p <- ifelse( data$DEL_maj_broad_chr_17p == 1 | data$DEL_maj_focal_TP53 == 1, 1, 0)
data$DEL_maj_call_chr_1p <- ifelse( data$DEL_maj_broad_chr_1p == 1 | data$DEL_maj_focal_CDKN2C == 1 | data$DEL_maj_focal_FAM46C == 1 | data$DEL_maj_focal_FAF1 == 1,  1, 0)
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
data$D_LAB_chem_albumin %>% summary

p1 <- data %>% ggplot(aes(MM_group, D_LAB_chem_albumin)) + 
  geom_violin() + 
  geom_boxplot(aes(fill=MM_group), alpha=0.5) + 
  stat_summary(geom="text", fun.y=quantile,
               aes(label=sprintf("%1.1f", ..y..), color=MM_group),
               position=position_nudge(x=0.48), size=3) +
  stat_compare_means(comparisons = my_comparisons)+ 
   stat_compare_means(label.y.npc = "bottom") + xlab("") 


#______ calcium _______ ns
data$D_LAB_chem_calcium %>% summary

p2 <- data %>% ggplot(aes(MM_group, D_LAB_chem_calcium)) + 
  geom_violin() + 
  geom_boxplot(aes(fill=MM_group), alpha=0.5) + 
  stat_summary(geom="text", fun.y=quantile,
               aes(label=sprintf("%1.1f", ..y..), color=MM_group),
               position=position_nudge(x=0.48), size=3) +
  stat_compare_means(comparisons = my_comparisons)+ 
   stat_compare_means(label.y.npc = "bottom") + xlab("") 


#______ creatinine _______ ns
data$D_LAB_chem_creatinine %>% summary

p3 <- data %>% ggplot(aes(MM_group, D_LAB_chem_creatinine)) + 
  geom_violin() + 
  geom_boxplot(aes(fill=MM_group), alpha=0.5) + 
  stat_summary(geom="text", fun.y=quantile,
               aes(label=sprintf("%1.1f", ..y..), color=MM_group),
               position=position_nudge(x=0.48), size=3) +
  stat_compare_means(comparisons = my_comparisons)+ 
   stat_compare_means(label.y.npc = "bottom") + xlab("") 


#______ hemoglobin _______ ***
data$D_LAB_cbc_hemoglobin %>% summary

p4 <- data %>% ggplot(aes(MM_group, D_LAB_cbc_hemoglobin)) + 
  geom_violin() + 
  geom_boxplot(aes(fill=MM_group), alpha=0.5) + 
  stat_summary(geom="text", fun.y=quantile,
               aes(label=sprintf("%1.1f", ..y..), color=MM_group),
               position=position_nudge(x=0.48), size=3) +
  stat_compare_means(comparisons = my_comparisons)+ 
   stat_compare_means(label.y.npc = "bottom") + xlab("") 


#______ LDH _______ ns
data$D_LAB_chem_ldh %>% summary

p5 <- data %>% ggplot(aes(MM_group, D_LAB_chem_ldh)) + 
  geom_violin() + 
  geom_boxplot(aes(fill=MM_group), alpha=0.5) + 
  stat_summary(geom="text", fun.y=quantile,
               aes(label=sprintf("%1.1f", ..y..), color=MM_group),
               position=position_nudge(x=0.48), size=3) +
  stat_compare_means(comparisons = my_comparisons)+ 
   stat_compare_means(label.y.npc = "bottom") + xlab("") 



#______ M-protein _______ ns
data$D_LAB_serum_m_protein %>% summary

p6 <- data %>% ggplot(aes(MM_group, D_LAB_serum_m_protein)) + 
  geom_violin() + 
  geom_boxplot(aes(fill=MM_group), alpha=0.5) + 
  stat_summary(geom="text", fun.y=quantile,
               aes(label=sprintf("%1.1f", ..y..), color=MM_group),
               position=position_nudge(x=0.48), size=3) +
  stat_compare_means(comparisons = my_comparisons)+ 
   stat_compare_means(label.y.npc = "bottom") + xlab("") 



#______ bun (blood urea nitrogen) _______ ns
data$D_LAB_chem_bun %>% summary

p7 <- data %>% ggplot(aes(MM_group, D_LAB_chem_bun)) + 
  geom_violin() + 
  geom_boxplot(aes(fill=MM_group), alpha=0.5) + 
  stat_summary(geom="text", fun.y=quantile,
               aes(label=sprintf("%1.1f", ..y..), color=MM_group),
               position=position_nudge(x=0.48), size=3) +
  stat_compare_means(comparisons = my_comparisons)+ 
   stat_compare_means(label.y.npc = "bottom") + xlab("") 



#______ abs neut _______ ns
data$D_LAB_cbc_abs_neut %>% summary

p8 <- data %>% ggplot(aes(MM_group, D_LAB_cbc_abs_neut)) + 
  geom_violin() + 
  geom_boxplot(aes(fill=MM_group), alpha=0.5) + 
  stat_summary(geom="text", fun.y=quantile,
               aes(label=sprintf("%1.1f", ..y..), color=MM_group),
               position=position_nudge(x=0.48), size=3) +
  stat_compare_means(comparisons = my_comparisons)+ 
   stat_compare_means(label.y.npc = "bottom") + xlab("") 



#______ glucose _______ ns
data$D_LAB_chem_glucose %>% summary

p9 <- data %>% ggplot(aes(MM_group, D_LAB_chem_glucose)) + 
  geom_violin() + 
  geom_boxplot(aes(fill=MM_group), alpha=0.5) + 
  stat_summary(geom="text", fun.y=quantile,
               aes(label=sprintf("%1.1f", ..y..), color=MM_group),
               position=position_nudge(x=0.48), size=3) +
  stat_compare_means(comparisons = my_comparisons)+ 
   stat_compare_means(label.y.npc = "bottom") + xlab("") 



#______ kappa _______ ns
data$D_LAB_serum_kappa %>% summary

p10 <- data %>% ggplot(aes(MM_group, D_LAB_serum_kappa)) + 
  geom_violin() + 
  geom_boxplot(aes(fill=MM_group), alpha=0.5) + 
  stat_summary(geom="text", fun.y=quantile,
               aes(label=sprintf("%1.1f", ..y..), color=MM_group),
               position=position_nudge(x=0.48), size=3) +
  stat_compare_means(comparisons = my_comparisons)+ 
   stat_compare_means(label.y.npc = "bottom") + xlab("") 




#______ lambda _______ ns
data$D_LAB_serum_lambda %>% summary

p11 <- data %>% ggplot(aes(MM_group, D_LAB_serum_lambda)) + 
  geom_violin() + 
  geom_boxplot(aes(fill=MM_group), alpha=0.5) + 
  stat_summary(geom="text", fun.y=quantile,
               aes(label=sprintf("%1.1f", ..y..), color=MM_group),
               position=position_nudge(x=0.48), size=3) +
  stat_compare_means(comparisons = my_comparisons)+ 
   stat_compare_means(label.y.npc = "bottom") + xlab("") 


#______ platlet _______ ***
data$D_LAB_cbc_platelet %>% summary

p12 <- data %>% ggplot(aes(MM_group, D_LAB_cbc_platelet)) + 
  geom_violin() + 
  geom_boxplot(aes(fill=MM_group), alpha=0.5) + 
  stat_summary(geom="text", fun.y=quantile,
               aes(label=sprintf("%1.1f", ..y..), color=MM_group),
               position=position_nudge(x=0.48), size=3) +
  stat_compare_means(comparisons = my_comparisons)+ 
   stat_compare_means(label.y.npc = "bottom") + xlab("") 



#______ totprot _______ ns
data$D_LAB_chem_totprot %>% summary

p13 <- data %>% ggplot(aes(MM_group, D_LAB_chem_totprot)) + 
  geom_violin() + 
  geom_boxplot(aes(fill=MM_group), alpha=0.5) + 
  stat_summary(geom="text", fun.y=quantile,
               aes(label=sprintf("%1.1f", ..y..), color=MM_group),
               position=position_nudge(x=0.48), size=3) +
  stat_compare_means(comparisons = my_comparisons)+ 
   stat_compare_means(label.y.npc = "bottom") + xlab("") 



#______ wbc _______ ns
data$D_LAB_cbc_wbc %>% summary

p14 <- data %>% ggplot(aes(MM_group, D_LAB_cbc_wbc)) + 
  geom_violin() + 
  geom_boxplot(aes(fill=MM_group), alpha=0.5) + 
  stat_summary(geom="text", fun.y=quantile,
               aes(label=sprintf("%1.1f", ..y..), color=MM_group),
               position=position_nudge(x=0.48), size=3) +
  stat_compare_means(comparisons = my_comparisons)+ 
   stat_compare_means(label.y.npc = "bottom") + xlab("") 



#______ beta2_microglobulin _______ .
data$D_LAB_serum_beta2_microglobulin %>% summary

p15 <- data %>% ggplot(aes(MM_group, D_LAB_serum_beta2_microglobulin)) + 
  geom_violin() + 
  geom_boxplot(aes(fill=MM_group), alpha=0.5) + 
  stat_summary(geom="text", fun.y=quantile,
               aes(label=sprintf("%1.1f", ..y..), color=MM_group),
               position=position_nudge(x=0.48), size=3) +
  stat_compare_means(comparisons = my_comparisons)+ 
   stat_compare_means(label.y.npc = "bottom") + xlab("") 




#______ c_reactive_protein _______ ns
data$D_LAB_serum_c_reactive_protein %>% summary

p16 <- data %>% ggplot(aes(MM_group, D_LAB_serum_c_reactive_protein)) + 
  geom_violin() + 
  geom_boxplot(aes(fill=MM_group), alpha=0.5) + 
  stat_summary(geom="text", fun.y=quantile,
               aes(label=sprintf("%1.1f", ..y..), color=MM_group),
               position=position_nudge(x=0.48), size=3) +
  stat_compare_means(comparisons = my_comparisons)+ 
  stat_compare_means(label.y.npc = "bottom") + xlab("")




Panel1 <- ggarrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,
          labels = c(paste0("a",1:16)),
          ncol = 3, nrow = 6)

Panel1 

ggsave("panel1.png", Panel1, path = "C:/Users/mm_gr/Desktop/", dpi = 300, width = 14, height = 25)


####################### panel 2 ##############################################


#______ IgA _______ ns
data$D_LAB_serum_iga %>% summary

q1 <- data %>% ggplot(aes(MM_group, D_LAB_serum_iga)) + 
  geom_violin() + 
  geom_boxplot(aes(fill=MM_group), alpha=0.5) + 
  stat_summary(geom="text", fun.y=quantile,
               aes(label=sprintf("%1.1f", ..y..), color=MM_group),
               position=position_nudge(x=0.48), size=3) +
  stat_compare_means(comparisons = my_comparisons)+ 
   stat_compare_means(label.y.npc = "bottom") + xlab("") 

#______ IgG _______ ns
data$D_LAB_serum_igg %>% summary

q2 <- data %>% ggplot(aes(MM_group, D_LAB_serum_igg)) + 
  geom_violin() + 
  geom_boxplot(aes(fill=MM_group), alpha=0.5) + 
  stat_summary(geom="text", fun.y=quantile,
               aes(label=sprintf("%1.1f", ..y..), color=MM_group),
               position=position_nudge(x=0.48), size=3) +
  stat_compare_means(comparisons = my_comparisons)+ 
   stat_compare_means(label.y.npc = "bottom") + xlab("") 

#______ IgE _______ ns
data$D_LAB_serum_ige %>% summary

q3 <- data %>% ggplot(aes(MM_group, D_LAB_serum_ige)) + 
  geom_violin() + 
  geom_boxplot(aes(fill=MM_group), alpha=0.5) + 
  stat_summary(geom="text", fun.y=quantile,
               aes(label=sprintf("%1.1f", ..y..), color=MM_group),
               position=position_nudge(x=0.48), size=3) +
  stat_compare_means(comparisons = my_comparisons)+ 
   stat_compare_means(label.y.npc = "bottom") + xlab("") 

#______ IgM _______ *
data$D_LAB_serum_igm %>% summary

q4 <- data %>% ggplot(aes(MM_group, D_LAB_serum_igm)) + 
  geom_violin() + 
  geom_boxplot(aes(fill=MM_group), alpha=0.5) + 
  stat_summary(geom="text", fun.y=quantile,
               aes(label=sprintf("%1.1f", ..y..), color=MM_group),
               position=position_nudge(x=0.48), size=3) +
  stat_compare_means(comparisons = my_comparisons)+ 
   stat_compare_means(label.y.npc = "bottom") + xlab("") 

#______ IgD _______ ns
data$D_LAB_serum_igd %>% summary

q5 <- data %>% ggplot(aes(MM_group, D_LAB_serum_igd)) + 
  geom_violin() + 
  geom_boxplot(aes(fill=MM_group), alpha=0.5) + 
  stat_summary(geom="text", fun.y=quantile,
               aes(label=sprintf("%1.1f", ..y..), color=MM_group),
               position=position_nudge(x=0.48), size=3) +
  stat_compare_means(comparisons = my_comparisons)+ 
   stat_compare_means(label.y.npc = "bottom") + xlab("") 

#______ CMMC _______ ***
data$CMMC %>% summary

q6 <- data %>% ggplot(aes(MM_group, CMMC)) + 
  geom_violin() + 
  geom_boxplot(aes(fill=MM_group), alpha=0.5) + 
  stat_summary(geom="text", fun.y=quantile,
               aes(label=sprintf("%1.1f", ..y..), color=MM_group),
               position=position_nudge(x=0.48), size=3) +
  stat_compare_means(comparisons = my_comparisons)+ 
  stat_compare_means(label.y = 70000) 




Panel2 <- ggarrange(q1,q2,q3,q4,q5,q6,
                    labels = c(paste0("b",1:6)),
                    ncol = 3, nrow = 6)

Panel2 

ggsave("panel2.png", Panel2, path = "C:/Users/mm_gr/Desktop/", dpi = 300, width = 14, height = 25)



#################### panel 3 ########################################



#______ CD38 ________ ns
data$D_IM_CD38_PC_PERCENT %>% summary

r1 <- data %>% ggplot(aes(MM_group, D_IM_CD38_PC_PERCENT)) + 
  geom_violin() + 
  geom_boxplot(aes(fill=MM_group), alpha=0.5) + 
  stat_summary(geom="text", fun.y=quantile,
               aes(label=sprintf("%1.1f", ..y..), color=MM_group),
               position=position_nudge(x=0.48), size=3) +
  stat_compare_means(comparisons = my_comparisons)+ 
   stat_compare_means(label.y.npc = "bottom") + xlab("") 


#______ CD138 _______ ns
data$D_IM_CD138_PC_PERCENT %>% summary

r2 <- data %>% ggplot(aes(MM_group, D_IM_CD138_PC_PERCENT)) + 
  geom_violin() + 
  geom_boxplot(aes(fill=MM_group), alpha=0.5) + 
  stat_summary(geom="text", fun.y=quantile,
               aes(label=sprintf("%1.1f", ..y..), color=MM_group),
               position=position_nudge(x=0.48), size=3) +
  stat_compare_means(comparisons = my_comparisons)+ 
   stat_compare_means(label.y.npc = "bottom") + xlab("") 

#______ CD56 _______ ns
data$D_IM_CD56_PC_PERCENT %>% summary

r3 <- data %>% ggplot(aes(MM_group, D_IM_CD56_PC_PERCENT)) + 
  geom_violin() + 
  geom_boxplot(aes(fill=MM_group), alpha=0.5) + 
  stat_summary(geom="text", fun.y=quantile,
               aes(label=sprintf("%1.1f", ..y..), color=MM_group),
               position=position_nudge(x=0.48), size=3) +
  stat_compare_means(comparisons = my_comparisons)+ 
   stat_compare_means(label.y.npc = "bottom") + xlab("") 

#______ CD13 _______ ns
data$D_IM_CD13_PC_PERCENT %>% summary

r4 <- data %>% ggplot(aes(MM_group, D_IM_CD13_PC_PERCENT)) + 
  geom_violin() + 
  geom_boxplot(aes(fill=MM_group), alpha=0.5) + 
  stat_summary(geom="text", fun.y=quantile,
               aes(label=sprintf("%1.1f", ..y..), color=MM_group),
               position=position_nudge(x=0.48), size=3) +
  stat_compare_means(comparisons = my_comparisons)+ 
   stat_compare_means(label.y.npc = "bottom") + xlab("") 

#______ CD19 _______ ns
data$D_IM_CD19_PC_PERCENT %>% summary

r5 <- data %>% ggplot(aes(MM_group, D_IM_CD19_PC_PERCENT)) + 
  geom_violin() + 
  geom_boxplot(aes(fill=MM_group), alpha=0.5) + 
  stat_summary(geom="text", fun.y=quantile,
               aes(label=sprintf("%1.1f", ..y..), color=MM_group),
               position=position_nudge(x=0.48), size=3) +
  stat_compare_means(comparisons = my_comparisons)+ 
   stat_compare_means(label.y.npc = "bottom") + xlab("") 

#______ CD20 _______ *
data$D_IM_CD20_PC_PERCENT %>% summary
table(data$D_IM_CD20_PC_PERCENT %>% is.na(), data$MM_group)

r6 <- data %>% ggplot(aes(MM_group, D_IM_CD20_PC_PERCENT)) + 
  geom_violin() + 
  geom_boxplot(aes(fill=MM_group), alpha=0.5) + 
  stat_summary(geom="text", fun.y=quantile,
               aes(label=sprintf("%1.1f", ..y..), color=MM_group),
               position=position_nudge(x=0.48), size=3) +
  stat_compare_means(comparisons = my_comparisons) +
                     # symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05,0.1, 1), symbols = c("****", "***", "**", "*", ".","ns")))+ 
   stat_compare_means(label.y.npc = "bottom") + xlab("") 

#______ CD27 _______ ns
data$D_IM_CD27_PC_PERCENT %>% summary
table(data$D_IM_CD27_PC_PERCENT %>% is.na(), data$MM_group)

r7 <- data %>% ggplot(aes(MM_group, D_IM_CD27_PC_PERCENT)) + 
  geom_violin() + 
  geom_boxplot(aes(fill=MM_group), alpha=0.5) + 
  stat_summary(geom="text", fun.y=quantile,
               aes(label=sprintf("%1.1f", ..y..), color=MM_group),
               position=position_nudge(x=0.48), size=3) +
  stat_compare_means(comparisons = my_comparisons)+ 
   stat_compare_means(label.y.npc = "bottom") + xlab("") 

#______ CD28 _______ ns
data$D_IM_CD28_PC_PERCENT %>% summary
table(data$D_IM_CD28_PC_PERCENT %>% is.na(), data$MM_group)

r8 <- data %>% ggplot(aes(MM_group, D_IM_CD28_PC_PERCENT)) + 
  geom_violin() + 
  geom_boxplot(aes(fill=MM_group), alpha=0.5) + 
  stat_summary(geom="text", fun.y=quantile,
               aes(label=sprintf("%1.1f", ..y..), color=MM_group),
               position=position_nudge(x=0.48), size=3) +
  stat_compare_means(comparisons = my_comparisons)+ 
   stat_compare_means(label.y.npc = "bottom") + xlab("") 


#______ CD33 _______ ns
data$D_IM_CD33_PC_PERCENT %>% summary
table(data$D_IM_CD33_PC_PERCENT %>% is.na(), data$MM_group)

r9 <- data %>% ggplot(aes(MM_group, D_IM_CD33_PC_PERCENT)) + 
  geom_violin() + 
  geom_boxplot(aes(fill=MM_group), alpha=0.5) + 
  stat_summary(geom="text", fun.y=quantile,
               aes(label=sprintf("%1.1f", ..y..), color=MM_group),
               position=position_nudge(x=0.48), size=3) +
  stat_compare_means(comparisons = my_comparisons)+ 
   stat_compare_means(label.y.npc = "bottom") + xlab("") 

#______ CD52 _______ .
data$D_IM_CD52_PC_PERCENT %>% summary
table(data$D_IM_CD52_PC_PERCENT %>% is.na(), data$MM_group)

r10 <- data %>% ggplot(aes(MM_group, D_IM_CD52_PC_PERCENT)) + 
  geom_violin() + 
  geom_boxplot(aes(fill=MM_group), alpha=0.5) + 
  stat_summary(geom="text", fun.y=quantile,
               aes(label=sprintf("%1.1f", ..y..), color=MM_group),
               position=position_nudge(x=0.48), size=3) +
  stat_compare_means(comparisons = my_comparisons)+ 
   stat_compare_means(label.y.npc = "bottom") + xlab("") 


#______ CD81 _______ nv
data$D_IM_CD81_PC_PERCENT %>% summary
table(data$D_IM_CD81_PC_PERCENT %>% is.na(), data$MM_group)

r11 <- data %>% ggplot(aes(MM_group, D_IM_CD81_PC_PERCENT)) + 
  geom_violin() + 
  geom_boxplot(aes(fill=MM_group), alpha=0.5) + 
  stat_summary(geom="text", fun.y=quantile,
               aes(label=sprintf("%1.1f", ..y..), color=MM_group),
               position=position_nudge(x=0.48), size=3) +
  stat_compare_means(comparisons = my_comparisons)+ 
   stat_compare_means(label.y.npc = "bottom") + xlab("") 


#______ CD117 _______ ***
data$D_IM_CD117_PC_PERCENT %>% summary
table(data$D_IM_CD117_PC_PERCENT %>% is.na(), data$MM_group)

r12 <- data %>% ggplot(aes(MM_group, D_IM_CD117_PC_PERCENT)) + 
  geom_violin() + 
  geom_boxplot(aes(fill=MM_group), alpha=0.5) + 
  stat_summary(geom="text", fun.y=quantile,
               aes(label=sprintf("%1.1f", ..y..), color=MM_group),
               position=position_nudge(x=0.48), size=3) +
  stat_compare_means(comparisons = my_comparisons)+ 
   stat_compare_means(label.y.npc = "bottom") + xlab("") 


#______ CD45 _______ ns
data$D_IM_CD45_PC_PERCENT %>% summary
table(data$D_IM_CD45_PC_PERCENT %>% is.na(), data$MM_group)

r13 <- data %>% ggplot(aes(MM_group, D_IM_CD45_PC_PERCENT)) + 
  geom_violin() + 
  geom_boxplot(aes(fill=MM_group), alpha=0.5) + 
  stat_summary(geom="text", fun.y=quantile,
               aes(label=sprintf("%1.1f", ..y..), color=MM_group),
               position=position_nudge(x=0.48), size=3) +
  stat_compare_means(comparisons = my_comparisons)+ 
  stat_compare_means(label.y.npc = "bottom") + xlab("") 


#______ CD319 _______ ns
data$D_IM_CD319_PC_PERCENT %>% summary
table(data$D_IM_CD319_PC_PERCENT %>% is.na(), data$MM_group)

r14 <- data %>% ggplot(aes(MM_group, D_IM_CD319_PC_PERCENT)) + 
  geom_violin() + 
  geom_boxplot(aes(fill=MM_group), alpha=0.5) + 
  stat_summary(geom="text", fun.y=quantile,
               aes(label=sprintf("%1.1f", ..y..), color=MM_group),
               position=position_nudge(x=0.48), size=3) +
  stat_compare_means(comparisons = my_comparisons)+ 
  stat_compare_means(label.y.npc = "bottom") + xlab("") 



#______ PERCENT_PC_IN_BM _______ ns
data$D_IM_MORPHOLOGY_PERCENT_PC_IN_BM %>% summary

r15 <- data %>% ggplot(aes(MM_group, D_IM_MORPHOLOGY_PERCENT_PC_IN_BM)) + 
  geom_violin() + 
  geom_boxplot(aes(fill=MM_group), alpha=0.5) + 
  stat_summary(geom="text", fun.y=quantile,
               aes(label=sprintf("%1.1f", ..y..), color=MM_group),
               position=position_nudge(x=0.48), size=3) +
  stat_compare_means(comparisons = my_comparisons)+ 
  stat_compare_means(label.y.npc = "bottom") + xlab("") 


#______ PERCENT_PC_IN_PB _______ ***
data$D_IM_MORPHOLOGY_PERCENT_PC_IN_PB %>% summary

r16 <- data %>% ggplot(aes(MM_group, D_IM_MORPHOLOGY_PERCENT_PC_IN_PB)) + 
  geom_violin() + 
  geom_boxplot(aes(fill=MM_group), alpha=0.5) + 
  stat_summary(geom="text", fun.y=quantile,
               aes(label=sprintf("%1.1f", ..y..), color=MM_group),
               position=position_nudge(x=0.48), size=3) +
  stat_compare_means(comparisons = my_comparisons)+ 
  stat_compare_means(label.y.npc = "bottom") + xlab("") 




Panel3 <- ggarrange(r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,r12,r13,r14,r15,r16,
          labels = c(paste0("c",1:16)),
          ncol = 3, nrow = 6)

Panel3

ggsave("panel3.png", Panel3, path = "C:/Users/mm_gr/Desktop/", dpi = 300, width = 14, height = 25)
