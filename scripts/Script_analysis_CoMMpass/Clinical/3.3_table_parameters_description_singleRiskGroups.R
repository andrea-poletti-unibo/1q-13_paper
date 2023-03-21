

############################ DEF TABLE in PAPER #############################

library(data.table)
library(tidyverse)
library(kableExtra)


tab_num_r1<- readxl::read_xlsx("C:/Users/andre/Alma Mater Studiorum Università di Bologna/PROJECT 1q & 13 - Documenti/Supplementary/Tables/clinical_variables_in_MMrisk1_2_3_CoMMpass.xlsx", sheet = 1)
tab_cat_r1<- readxl::read_xlsx("C:/Users/andre/Alma Mater Studiorum Università di Bologna/PROJECT 1q & 13 - Documenti/Supplementary/Tables/clinical_variables_in_MMrisk1_2_3_CoMMpass.xlsx", sheet = 2)

tab_num_r2<- readxl::read_xlsx("C:/Users/andre/Alma Mater Studiorum Università di Bologna/PROJECT 1q & 13 - Documenti/Supplementary/Tables/clinical_variables_in_MMrisk1_2_3_CoMMpass.xlsx", sheet = 3)
tab_cat_r2<- readxl::read_xlsx("C:/Users/andre/Alma Mater Studiorum Università di Bologna/PROJECT 1q & 13 - Documenti/Supplementary/Tables/clinical_variables_in_MMrisk1_2_3_CoMMpass.xlsx", sheet = 4)

tab_num_r3<- readxl::read_xlsx("C:/Users/andre/Alma Mater Studiorum Università di Bologna/PROJECT 1q & 13 - Documenti/Supplementary/Tables/clinical_variables_in_MMrisk1_2_3_CoMMpass.xlsx", sheet = 5)
tab_cat_r3<- readxl::read_xlsx("C:/Users/andre/Alma Mater Studiorum Università di Bologna/PROJECT 1q & 13 - Documenti/Supplementary/Tables/clinical_variables_in_MMrisk1_2_3_CoMMpass.xlsx", sheet = 6)

tab_cat_r1$variable

# hyperdiploidy
# 1p36 del (FISH)
# 17p del (FISH)
# t(4;14) (FISH)
# t(14;16) (FISH)
# CN gain 11p
# CN gain 11q
# CN loss 1q
# CN gain 18p
# CN gain 18q
# MYC CN gain (focal)
# cdkn2c CN loss (focal)
# CYLD CN loss (focal)
# FAF1 CN loss (focal)
# FAM46C CN loss (focal)
# TP53 CN loss (focal)
# TRAF3 CN loss (focal)
# YAP1/BIRC2/3 CN loss (focal)
# NOTCH2 CN gain (focal)

VARS <- c("HyperDiploidy", 
          "SeqWGS_WHSC1_CALL", 
          "SeqWGS_CCND1_CALL",
          "SeqWGS_CCND3_CALL",
          "SeqWGS_MAF_CALL", 
          "SeqWGS_MAFB_CALL",
          "SeqWGS_MYC_CALL",
          "SeqWGS_CCND2_CALL",
          "SeqWGS_MAFA_CALL",
          "AMP_maj_broad_chr_11p",
          "AMP_maj_broad_chr_11q",
          "DEL_maj_broad_chr_1p",
          "AMP_maj_broad_chr_18p",
          "AMP_maj_broad_chr_18q",
          "AMP_maj_focal_MYC",
          "DEL_maj_focal_TP53",
          "DEL_maj_focal_CDKN2C",
          "DEL_maj_focal_FAM46C",
          "DEL_maj_focal_CYLD",
          "DEL_maj_focal_TRAF3")

# tab_cat_r1_VARS <- tab_cat_r1 %>% filter(type== 1, variable %in% VARS) %>% select(variable, `1q&13+ perc`=`1q&13+_perc`, "1q&13+ other perc"= other_perc, "1q&13+ p.value"= p.value)
# tab_cat_r2_VARS <- tab_cat_r2 %>% filter(type== 1, variable %in% VARS) %>% select( `1q/13 perc`=`1q/13_perc`, "1q/13 other perc"= other_perc, "1q/13 p.value"= p.value)
# tab_cat_r3_VARS <- tab_cat_r3 %>% filter(type== 1, variable %in% VARS) %>% select( `1q&13- perc`=`1q&13-_perc`, "1q&13- other perc"= other_perc, "1q&13- p.value"= p.value)


tab_cat_r1_VARS <- tab_cat_r1 %>% filter(type== 1, variable %in% VARS) %>% select(variable, `1q&13+`, other, `1q&13+ perc`=`1q&13+_perc`, "other perc"= other_perc, "p.value"= p.value)
tab_cat_r1_VARS <- tab_cat_r1_VARS %>% mutate("code" = case_when(p.value>0.1 ~ "n.s.",
                                                                 p.value>0.05 ~ ".", 
                                                                 p.value>0.01 ~ "*",
                                                                 p.value>0.001 ~ "**",
                                                                 p.value<= 0.001 ~ "***")) %>% 
  mutate(p.value = p.value %>% round(3)) %>% 
  mutate(variable=variable %>% str_remove("maj-broad_|maj-focal_")) %>% 
  mutate(variable=variable %>% str_replace("SeqWGS_(.*)_CALL","traslocation \\1")) %>% 
  mutate(variable=variable %>% str_replace_all("_"," ")) %>% 
  mutate(`1q&13+ perc` = `1q&13+ perc` %>% `*`(100) %>% paste0(.,"%")) %>% 
  mutate(`other perc` = `other perc` %>% `*`(100) %>% paste0(.,"%")) %>% 
  mutate(`1q&13+` = paste0(`1q&13+`," (", `1q&13+ perc`,")")) %>% 
  mutate(`other` = paste0(`other`," (", `other perc`,")")) %>% 
  select(-contains("perc"))


tab_cat_r2_VARS <- tab_cat_r2 %>% filter(type== 1, variable %in% VARS) %>% select(`1q/13`, other, `1q/13 perc`=`1q/13_perc`, "other perc"= other_perc, "p.value"= p.value)
tab_cat_r2_VARS <- tab_cat_r2_VARS %>% mutate("code" = case_when(p.value>0.1 ~ "n.s.",
                                                                 p.value>0.05 ~ ".", 
                                                                 p.value>0.01 ~ "*",
                                                                 p.value>0.001 ~ "**",
                                                                 p.value<= 0.001 ~ "***")) %>% 
  mutate(p.value = p.value %>% round(3)) %>% 
  mutate(`1q/13 perc` = `1q/13 perc` %>% `*`(100) %>% paste0(.,"%")) %>% 
  mutate(`other perc` = `other perc` %>% `*`(100) %>% paste0(.,"%")) %>% 
  mutate(`1q/13` = paste0(`1q/13`," (", `1q/13 perc`,")")) %>% 
  mutate(`other` = paste0(`other`," (", `other perc`,")")) %>% 
  select(-contains("perc"))


tab_cat_r3_VARS <- tab_cat_r3 %>% filter(type== 1, variable %in% VARS) %>% select( `1q&13-`, other, `1q&13- perc`=`1q&13-_perc`, "other perc"= other_perc, "p.value"= p.value)
tab_cat_r3_VARS <- tab_cat_r3_VARS %>% mutate("code" = case_when(p.value>0.1 ~ "n.s.",
                                                                 p.value>0.05 ~ ".", 
                                                                 p.value>0.01 ~ "*",
                                                                 p.value>0.001 ~ "**",
                                                                 p.value<= 0.001 ~ "***")) %>% 
  mutate(p.value = p.value %>% round(3)) %>% 
  mutate(`1q&13- perc` = `1q&13- perc` %>% `*`(100) %>% paste0(.,"%")) %>% 
  mutate(`other perc` = `other perc` %>% `*`(100) %>% paste0(.,"%")) %>% 
  mutate(`1q&13-` = paste0(`1q&13-`," (", `1q&13- perc`,")")) %>% 
  mutate(`other` = paste0(`other`," (", `other perc`,")")) %>% 
  select(-contains("perc"))


all <- cbind(tab_cat_r1_VARS, tab_cat_r2_VARS, tab_cat_r3_VARS)


all %>% kbl(caption = "1q&13 groups desctiption table - CoMMpass dataset") %>%
  kable_classic(full_width = F, html_font = "Cambria") %>% save_kable("table_risk123_description_CoMMpass.pdf",vwidth=500, vheight=300)


