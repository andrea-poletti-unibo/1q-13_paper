
library(tidyverse)
library(data.table)

# data <- fread("C:/Users/mm_gr/Alma Mater Studiorum Università di Bologna/PROJECT 1q & 13 - Documenti/complete_database_1q_13_CoMMpass_240121.txt", na.strings=c("","NA"))


import <- fread("D:/analisi_in_corso/risk_1q_13/CoMMpass/complete_database_1q_13_CoMMpass_240121.txt") 


# add sex and clincal data
import_clinical <- fread("D:/CoMMpass_data/IA13a/Clinical_Flat_Files/CoMMpass_IA13_FlatFiles/MMRF_CoMMpass_IA13_PER_PATIENT.csv")

import_clinical$Study_Visit_iD <- paste0(import_clinical$PUBLIC_ID, "_1_BM")

import$Study_Visit_iD
import_clinical$Study_Visit_iD

import2 <- left_join(import, import_clinical,  by = "Study_Visit_iD")

import <- import2

data <- import


data$PUBLIC_ID <- data$Study_Visit_iD

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

#_____ create tables _______

data %>% names

data <- data %>% as.data.frame() 

ID <- "PUBLIC_ID"
groupvar <- "MM_group"

groups <- data[,groupvar] %>% unique() %>% as.character()


################## numerical variables ###################

data$PERFORMANCE <- parse_number(data$ECOG_PERFORMANCEST)
data$MIN_LYTICLESION <- parse_number(data$BA_OFLYTICLESION)


data_var_num <- data %>% select(where(function(x) is.numeric(x) && sum(is.na(x))<nrow(data)*0.8 && length(unique(x))>2), 
                                !!as.symbol(groupvar), 
                                !!as.symbol(ID))



res <- NULL

i=data_var_num %>% names %>% .[5]

for(i in data_var_num %>% names){
  print(i)
  
  d <- data_var_num %>% 
    group_by(!!as.symbol(groupvar)) %>% 
    summarise(mean=mean(!!as.symbol(i), na.rm = T), 
              NAs=sum(is.na(!!as.symbol(i))), 
              n=sum(!is.na(!!as.symbol(i))))
  
  means <- d$mean
  NAs <- d$NAs
  n <- d$n
  
  try(krustal.pval <- kruskal.test(data_var_num[,i]~data_var_num[,groupvar])$p.value %>% round(4))
  
  name_means_p <- c(i, krustal.pval, means, n, NAs)

  res <- rbind(res, name_means_p) 
}

colnames(res) <- c("variable", 
                   "p.val",
                   paste0("mean_",d$MM_group %>% as.character()), 
                   paste0("n_",d$MM_group %>% as.character()),
                   paste0("NAs_",d$MM_group %>% as.character())
                   )

res_num_df <- res %>% as.data.frame()


xlsx::write.xlsx(res_num_df, 
                 "C:/Users/andre/Alma Mater Studiorum Università di Bologna/PROJECT 1q & 13 - Documenti/Supplementary/clinical_numeric_variables_in_MMgroups_Commpass.xlsx",
                 row.names = F)



################ categorical variables #######################

library(fastDummies)

data_var_cat <- data %>% select(where(function(x) !is.na(x) && length(unique(x))<15 && length(unique(x))>1), !!groupvar)

data_var_cat %>% names
data_var_cat_def <- data_var_cat %>% select(# !starts_with("D_TRI") & 
                                            !ends_with("DETECTED") &
                                            # !starts_with("HD") & 
                                            # !starts_with("DEL") & 
                                            !starts_with("D_IM"),
                                            # !starts_with("D_QOL"),
                                            !!groupvar
                                          )
data_var_cat_def %>% names

data_var_cat_def_dummy <- data_var_cat_def %>% dummy_cols(ignore_na = T, remove_selected_columns = F)

iter <- data_var_cat_def_dummy %>% names
iter

i <- "SEX_F"

res <- NULL

for(i in iter) {
  
  print(i)
  
  try(pval <- with(data_var_cat_def_dummy, {fisher.test(get(groupvar), get(i) )})$p.value %>% round(4))
  
  try(df <- with(data_var_cat_def_dummy, {table(get(groupvar), get(i)) %>% as.data.frame.matrix() %>% t %>% as.data.frame() %>% rownames_to_column(var = "type")}) )
  
  try(name_df <- cbind(variable=i, df, p.value=pval))
  
  try(res <- rbind(res, name_df) )
  
}

res_cat_df <- res %>% as.data.frame()


data$MM_group %>% table
res_cat_df$`1q&13+_perc` <- round(res_cat_df$`1q&13+`/187,3)
res_cat_df$`1q/13_perc`  <- round(res_cat_df$`1q/13`/317,3)
res_cat_df$`1q&13-_perc` <- round(res_cat_df$`1q&13-`/336,3)


xlsx::write.xlsx(res_cat_df, 
                 "C:/Users/andre//Alma Mater Studiorum Università di Bologna/PROJECT 1q & 13 - Documenti/Supplementary/clinical_categoric_variables_in_MMgroups_CoMMpass.xlsx",
                 row.names = F)

