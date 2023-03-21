
library(data.table)
library(tidyverse)
library(RODBC)

import2 <- data.table::fread("C:/Users/mm_gr/Alma Mater Studiorum Università di Bologna/PROJECT 1q & 13 - Documenti/complete_database_1q_13_181019.txt")


db_clin <- RODBC::odbcConnectAccess2007("C:/Users/mm_gr/Alma Mater Studiorum Università di Bologna/MM group - Documenti/Databases/MM_db_clinical_records.accdb")
clin <- sqlFetch(db_clin, "Extraction_1q13_301120", dec=",", na.strings="nv" )


db_exp<- RODBC::odbcConnectAccess2007("C:/Users/mm_gr/Alma Mater Studiorum Università di Bologna/MM group - Documenti/Databases/MM_db_experiments.accdb")
calls <- sqlFetch(db_exp, "SNP_calls_1q_13_181019", dec=",")


import <- inner_join(clin, calls, by=c("SNP_ARRAY_ESORDIO"="SNP"))


#______ MM RISK CREATION _____
import$AMP_1q_genes_all <- with(import, 
                                ifelse(`AMP_maj-focal_ANP32E`==1 | `AMP_maj-focal_MCL1` ==1 | `AMP_maj-focal_CKS1B`==1, 1,0 ))


table(import$`AMP_maj-broad_chr_1q`, import$AMP_1q_genes_all)
table(import$`DEL_maj-broad_chr_13q`, import$`DEL_maj-focal_RB1`)


import$MMrisk_1q_all <- ifelse( import$`AMP_maj-broad_chr_1q` == 1 | import$AMP_1q_genes_all ==1, 1, 0)
import$MMrisk_13_all <- ifelse( import$`DEL_maj-broad_chr_13q` ==1 | import$`DEL_maj-focal_RB1` ==1, 1, 0)

import$MMrisk_CLASSIFIER_ALL <- 3  - import$MMrisk_1q_all - import$MMrisk_13_all

table(import$MMrisk_CLASSIFIER_ALL)

import$MMrisk_AB_ALL <- ifelse(import$MMrisk_CLASSIFIER_ALL == 2 & import$MMrisk_1q_all ==1,
                               "2_1q", 
                               ifelse(import$MMrisk_CLASSIFIER_ALL == 2 & import$MMrisk_13_all ==1,
                                      "2_13",
                                      import$MMrisk_CLASSIFIER_ALL)
)


############################################

library(tidyverse)
library(ggvenn)

# 4 way VENN PLOT

is.na(venndf$TRANSLOCATION) %>% sum

venndf <- import %>% select(MMrisk_1q_all, MMrisk_13_all, HyperDiploidy, TRANSLOCATION)

venndf <- apply(venndf, 2 , as.logical) %>% as.data.frame() %>% rownames_to_column()

venndf_complete <- venndf[complete.cases(venndf),]

# venndf$TRANSLOCATION[is.na(venndf$TRANSLOCATION)] <- F

venndf_complete <- venndf_complete %>% rename("Amp 1q"=MMrisk_1q_all, "Del 13q"=MMrisk_13_all, "IgH Translocation"=TRANSLOCATION)

venndf_complete %>% names

venndf_complete <- venndf_complete %>% select(`Amp 1q`,HyperDiploidy,`IgH Translocation`,`Del 13q`)

ggvenn(venndf_complete, 
       fill_color = c("#1C86EE","#458B00", "#BF3EFF","#FF3030"),
       fill_alpha = 0.5,
       set_name_size = 5) + 
  annotate(geom = "text", label="BO dataset", x=-1.64, y=1.6, size=5)

ggsave("C:/Users/mm_gr/Alma Mater Studiorum Università di Bologna/PROJECT 1q & 13 - Documenti/PAPER_FIGURES_OFFICIAL/Venn4Ways/Venn_diagram_BO.pdf", 
       width = 7, height = 7)



