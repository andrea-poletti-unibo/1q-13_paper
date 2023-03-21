library(tidyverse)
library(RODBC)

db_clin <- RODBC::odbcConnectAccess2007("C:/Users/andre/Alma Mater Studiorum Università di Bologna/MM group - Documenti/EVERGREEN_inter-projects_files/MM_db_clinical_records.accdb")
clin <- sqlFetch(db_clin, "Extraction_1q13_301120", dec=",", na.strings="nv" )


db_exp<- RODBC::odbcConnectAccess2007("C:/Users/andre/Alma Mater Studiorum Università di Bologna/MM group - Documenti/EVERGREEN_inter-projects_files/MM_db_experiments.accdb")
calls <- sqlFetch(db_exp, "SNP_calls_1q_13_181019", dec=",")


data <- inner_join(clin, calls, by=c("SNP_ARRAY_ESORDIO"="SNP"))


#____ Build CALLS for arms with broad+focal_____

#amps
data$`AMP_maj-call_chr_1q` <- ifelse( data$`AMP_maj-broad_chr_1q` == 1 | data$`AMP_maj-focal_ANP32E`== 1 | data$`AMP_maj-focal_CKS1B` ==1 | data$`AMP_maj-focal_MCL1` ==1, 1, 0)
data$`AMP_maj-call_chr_8q` <- ifelse( data$`AMP_maj-broad_chr_8q` == 1 | data$`AMP_maj-focal_MYC` == 1, 1, 0)

#dels
data$`DEL_maj-call_chr_17p` <- ifelse( data$`DEL_maj-broad_chr_17p` == 1 | data$`DEL_maj-focal_TP53` == 1, 1, 0)
data$`DEL_maj-call_chr_1p`<-  ifelse( data$`DEL_maj-broad_chr_1p`== 1 |  data$`DEL_maj-focal_CDKN2C` == 1 | data$`DEL_maj-focal_FAM46C` == 1 | data$`DEL_maj-focal_FAF1` == 1,  1, 0)
data$`DEL_maj-call_chr_13q` <- ifelse( data$`DEL_maj-broad_chr_13q` == 1 | data$`DEL_maj-focal_RB1` == 1, 1, 0)
data$`DEL_maj-call_chr_14q` <- ifelse( data$`DEL_maj-broad_chr_14q` == 1 | data$`DEL_maj-focal_TRAF3` == 1, 1, 0)
data$`DEL_maj-call_chr_16q` <- ifelse( data$`DEL_maj-broad_chr_16q` == 1 | data$`DEL_maj-focal_CYLD` == 1, 1, 0)

#_____ add group risk 1q & 13 _______
MM_risk <- ifelse( data$`AMP_maj-call_chr_1q` == 1 & data$`DEL_maj-call_chr_13q` ==1, 1,
                   ifelse(data$`AMP_maj-call_chr_1q` == 1 | data$`DEL_maj-call_chr_13q` ==1, 2, 3) ) 

data$MM_group_num <- MM_risk
data$MM_group <- MM_risk %>% factor(levels = c("1","2","3"), labels = c("1q&13+","1q/13","1q&13-"))

data$MM_risk_1 <- ifelse(MM_risk==1, 1,0)
data$MM_risk_2 <- ifelse(MM_risk==2, 1,0)
data$MM_risk_3 <- ifelse(MM_risk==3, 1,0)

#_____ create tables _______

data %>% names

data <- data %>% as.data.frame() 

ID <- "Project_1q_13_Classificator_UPN"
groupvar <- "MM_group"


data$PROTOCOL_group <- ifelse(data$PROTOCOL=="EMN02","EMN02",ifelse(data$PROTOCOL=="BO2005","BO2005","AMB"))
data$PROTOCOL_group %>% table

groupvar <- "PROTOCOL_group"

groups <- data[,groupvar] %>% unique() %>% as.character()


################## numerical variables ###################


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
  
  try(name_means_p <- c(i, krustal.pval, means, n, NAs))
  
  try(res <- rbind(res, name_means_p))
}

colnames(res) <- c("variable", 
                   "p.val",
                   paste0("mean_",d[,1] %>% unlist()), 
                   paste0("n_",d[,1] %>% unlist()),
                   paste0("NAs_",d[,1] %>% unlist()))

res_num_df <- res %>% as.data.frame()


xlsx::write.xlsx(res_num_df, 
                 "C:/Users/mm_gr//Alma Mater Studiorum Università di Bologna/PROJECT 1q & 13 - Documenti/Supplementary/clinical_numeric_variables_in_MMgroups_Bologna.xlsx",
                 row.names = F)


xlsx::write.xlsx(res_num_df, 
                 "C:/Users/andre/Alma Mater Studiorum Università di Bologna/PROJECT 1q & 13 - Documenti/Supplementary/clinical_numeric_variables_in_PROTOCOLS_Bologna.xlsx",
                 row.names = F)


################ categorical variables #######################

library(fastDummies)

data_var_cat <- data %>% select(where(function(x) !is.na(x) && length(unique(x))<15 && length(unique(x))>1), 
                                !!as.symbol(groupvar), 
                                !!as.symbol(ID))

data_var_cat %>% names
data_var_cat_def <- data_var_cat %>% select(PROTOCOL:SEX, 
                                            FISH_T_4_14:LIGHT_CHAIN, 
                                            TX_I_LINE_1_0, 
                                            NUMBERS_TX_FRONT_LINE, 
                                            MAINTENANCE_YES_NO:MAINTENANCE_TERAPHY, 
                                            CONSOLIDATION_YES_NO,
                                            !!groupvar)

data_var_cat_def_dummy <- data_var_cat_def %>% dummy_cols(ignore_na = T)



i <- "SEX_F"

res <- NULL

for(i in data_var_cat_def_dummy %>% names) {
  
  print(i)

try(pval <- with(data_var_cat_def_dummy, {fisher.test(get(groupvar), get(i) )})$p.value %>% round(4))

try(df <- with(data_var_cat_def_dummy, {table(get(groupvar), get(i)) %>% as.data.frame.matrix() %>% t %>% as.data.frame() %>% rownames_to_column(var = "type")}) )

try(name_df <- cbind(variable=i, df, p.value=pval))

try(res <- rbind(res, name_df) )

}

res_cat_df <- res %>% as.data.frame()

xlsx::write.xlsx(res_cat_df, 
                 "C:/Users/mm_gr//Alma Mater Studiorum Università di Bologna/PROJECT 1q & 13 - Documenti/Supplementary/clinical_categoric_variables_in_MMgroups_Bologna.xlsx",
                 row.names = F)


xlsx::write.xlsx(res_cat_df, 
                 "C:/Users/andre/Alma Mater Studiorum Università di Bologna/PROJECT 1q & 13 - Documenti/Supplementary/clinical_categoric_variables_in_PROTOCOLS_Bologna.xlsx",
                 row.names = F)
