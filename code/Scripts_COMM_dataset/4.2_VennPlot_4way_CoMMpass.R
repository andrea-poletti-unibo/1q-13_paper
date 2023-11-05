
library(data.table)
library(tidyverse)

import <-fread("results/complete_database_CoMMpass_1q_13.txt")

#______ MM RISK CREATION _____
import$AMP_1q_genes_all <- with(import, 
                                ifelse(`AMP_maj_focal_ANP32E`==1 | `AMP_maj_focal_MCL1` ==1 | `AMP_maj_focal_CKS1B`==1, 1,0 ))


import$MMrisk_1q_all <- ifelse( import$`AMP_maj_broad_chr_1q` == 1 | import$AMP_1q_genes_all ==1, 1, 0)
import$MMrisk_13_all <- ifelse( import$`DEL_maj_broad_chr_13q` ==1 | import$`DEL_maj_focal_RB1` ==1, 1, 0)

import$MMrisk_CLASSIFIER_ALL <- 3  - import$MMrisk_1q_all - import$MMrisk_13_all

table(import$MMrisk_CLASSIFIER_ALL)

import$MMrisk_AB_ALL <- ifelse(import$MMrisk_CLASSIFIER_ALL == 2 & import$MMrisk_1q_all ==1,
                               "2_1q", 
                               ifelse(import$MMrisk_CLASSIFIER_ALL == 2 & import$MMrisk_13_all ==1,
                                      "2_13",
                                      import$MMrisk_CLASSIFIER_ALL)
)

#_______ TRANSLOCATION CREATION _________

import$TRANSLOCATION <- ifelse(import$SeqWGS_WHSC1_CALL==1 |
                                 import$SeqWGS_CCND1_CALL==1 |
                                 import$SeqWGS_CCND3_CALL==1 |
                                 import$SeqWGS_MAF_CALL ==1 |
                                 import$SeqWGS_MAFB_CALL ==1 |
                                 import$SeqWGS_CCND2_CALL ==1 |
                                 import$SeqWGS_MAFA_CALL ==1, 1, 0)


############################################

library(tidyverse)
library(ggvenn)

# 4 way VENN PLOT

venndf <- import %>% select(MMrisk_1q_all, MMrisk_13_all, HyperDiploidy, TRANSLOCATION)

venndf <- apply(venndf, 2 , as.logical) %>% as.data.frame() %>% rownames_to_column()

venndf_complete <- venndf[complete.cases(venndf),]

# venndf$TRANSLOCATION[is.na(venndf$TRANSLOCATION)] <- F

venndf_complete <- venndf_complete %>% rename("Amp 1q"=MMrisk_1q_all, "Del 13q"=MMrisk_13_all, "IgH Translocation"=TRANSLOCATION)

venndf_complete %>% names

venndf_complete <- venndf_complete %>% select(`Amp 1q`,HyperDiploidy,`IgH Translocation`,`Del 13q`)

ggvenn(venndf_complete, 
       fill_color = c("#1C86EE","#458B00","#BF3EFF","#FF3030"),
       fill_alpha = 0.5,
       set_name_size = 5) + 
  annotate(geom = "text", label="CoMMpass dataset", x=-1.44, y=1.6, size=5)


dir.create("plots/Venn4Ways/", showWarnings = F, recursive = T)

ggsave("plots/Venn4Ways/Venn_diagram_CoMM.pdf", 
       width = 7, height = 7)


ggsave("plots/Venn4Ways/Venn_diagram_CoMM.svg", 
       width = 7, height = 7)
