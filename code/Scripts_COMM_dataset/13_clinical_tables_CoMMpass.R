library(tidyverse)
library(data.table)

data <- fread("results/complete_database_CoMMpass_1q_13.txt")

dir.create("results/Clinical_anlysis/CoMM_dataset/tables/", recursive = T, showWarnings = F)


data$PUBLIC_ID
data$Study_Visit_iD


#_____ add group risk 1q & 13 _______
MM_risk <- ifelse( data$AMP_maj_call_chr_1q == 1 & data$DEL_maj_call_chr_13q ==1, 1,
                   ifelse(data$AMP_maj_call_chr_1q == 1 | data$DEL_maj_call_chr_13q ==1, 2, 3) ) 

data$MM_group_num <- MM_risk
data$MM_group <- MM_risk %>% factor(levels = c("1","2","3"), labels = c("1q&13+","1q/13","1q&13-"))

data$MM_risk_1 <- ifelse(MM_risk==1, 1,0)
data$MM_risk_2 <- ifelse(MM_risk==2, 1,0)
data$MM_risk_3 <- ifelse(MM_risk==3, 1,0)


data$PERFORMANCE <- parse_number(data$ECOG_PERFORMANCEST)
data$MIN_LYTICLESION <- parse_number(data$BA_OFLYTICLESION)



#_____ create tables _______

data %>% names

data <- data %>% as.data.frame() 

ID <- "PUBLIC_ID"




############################ MM_group variable #################################


groupvar <- "MM_group"

groups <- data[,groupvar] %>% unique() %>% as.character()


#============ numerical ============


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
                 "results/Clinical_anlysis/CoMM_dataset/tables/clinical_numeric_variables_in_MMgroups_CoMMpass.xlsx",
                 row.names = F)




#============ categorical ============


library(fastDummies)

data_var_cat <- data %>% select(where(function(x) !is.na(x) && length(unique(x))<15 && length(unique(x))>1), 
                                !!as.symbol(groupvar), 
                                !!as.symbol(ID))

data_var_cat %>% names
data_var_cat_def <- data_var_cat %>% select(# !starts_with("D_TRI") & 
  !ends_with("DETECTED") &
    # !starts_with("HD") & 
    # !starts_with("DEL") & 
    !starts_with("D_IM"),
  # !starts_with("D_QOL"),
  !!groupvar,
) %>% select(- PUBLIC_ID)

data_var_cat_def_dummy <- data_var_cat_def %>% dummy_cols(ignore_na = T)

i <- "SEX_F"

res <- NULL

for(i in data_var_cat_def_dummy %>% names) {
  
  print(i)
  
  try(pval <- with(data_var_cat_def_dummy, {fisher.test(get(groupvar), get(i) )})$p.value %>% round(4))
  
  try(GM <- with(data_var_cat_def_dummy, {gmodels::CrossTable(get(groupvar), get(i)) }) ) 
  
  try(df <- GM$t %>% as.data.frame.matrix() %>% t %>% as.data.frame() %>% rownames_to_column(var = "type")) 
  
  try(freq <- GM$prop.row %>% as.data.frame.matrix() %>% t %>% as.data.frame() %>% set_names(paste0(names(.), "_freq")) ) 
  
  try(tot <- GM$t %>% as.data.frame.matrix() %>% t %>% as.data.frame() %>% rowSums()) 
  
  try(all <- tot %>% sum) 
  
  try(tot_freq <- round(tot/all, 4)) 
  
  try(name_df <- cbind(variable=i, df, p.value=pval, freq, tot, tot_freq))
  
  try(res <- rbind(res, name_df) )
  
}

res_cat_df <- res %>% as.data.frame()

xlsx::write.xlsx(res_cat_df, 
                 "results/Clinical_anlysis/CoMM_dataset/tables/clinical_categoric_variables_in_MMgroups_CoMMpass.xlsx",
                 row.names = F)





############################ Risk 1 "1q&13+" variable #################################


data$MM_risk_1 <- data$MM_risk_1 %>% recode("1"= "1q&13+", "0"="other")

groupvar <- "MM_risk_1"
data$MM_risk_1

groups <- data[,groupvar] %>% unique() %>% as.character()



#============ numerical ============


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
                   paste0("mean_",d[,1] %>% unlist()), 
                   paste0("n_",d[,1] %>% unlist()),
                   paste0("NAs_",d[,1] %>% unlist()))


res_num_df <- res %>% as.data.frame()

xlsx::write.xlsx(res_num_df, 
                 "results/Clinical_anlysis/CoMM_dataset/tables/clinical_numeric_variables_in_group1_CoMMpass.xlsx",
                 row.names = F)




#============ categorical ============


library(fastDummies)

data_var_cat <- data %>% select(where(function(x) !is.na(x) && length(unique(x))<15 && length(unique(x))>1), 
                                !!as.symbol(groupvar), 
                                !!as.symbol(ID))

data_var_cat %>% names
data_var_cat_def <- data_var_cat %>% select(# !starts_with("D_TRI") & 
  !ends_with("DETECTED") &
    # !starts_with("HD") & 
    # !starts_with("DEL") & 
    !starts_with("D_IM"),
  # !starts_with("D_QOL"),
  !!groupvar,
) %>% select(- PUBLIC_ID)

data_var_cat_def_dummy <- data_var_cat_def %>% dummy_cols(ignore_na = T)

i <- "SEX_F"

res <- NULL

for(i in data_var_cat_def_dummy %>% names) {
  
  print(i)
  
  try(pval <- with(data_var_cat_def_dummy, {fisher.test(get(groupvar), get(i) )})$p.value %>% round(4))
  
  try(GM <- with(data_var_cat_def_dummy, {gmodels::CrossTable(get(groupvar), get(i)) }) ) 
  
  try(df <- GM$t %>% as.data.frame.matrix() %>% t %>% as.data.frame() %>% rownames_to_column(var = "type")) 
  
  try(freq <- GM$prop.row %>% as.data.frame.matrix() %>% t %>% as.data.frame() %>% set_names(paste0(names(.), "_freq")) ) 
  
  try(tot <- GM$t %>% as.data.frame.matrix() %>% t %>% as.data.frame() %>% rowSums()) 
  
  try(all <- tot %>% sum) 
  
  try(tot_freq <- round(tot/all, 4)) 
  
  try(name_df <- cbind(variable=i, df, p.value=pval, freq, tot, tot_freq))
  
  try(res <- rbind(res, name_df) )
  
}

res_cat_df <- res %>% as.data.frame()

xlsx::write.xlsx(res_cat_df, 
                 "results/Clinical_anlysis/CoMM_dataset/tables/clinical_categoric_variables_in_group1_CoMMpass.xlsx",
                 row.names = F)





############################ risk 3 "1q&13-" variable #################################


data$MM_risk_3 <- data$MM_risk_3 %>% recode("1"= "1q&13-", "0"="other")

groupvar <- "MM_risk_3"
data$MM_risk_3

groups <- data[,groupvar] %>% unique() %>% as.character()


#============ numerical ============


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
                   paste0("mean_",d[,1] %>% unlist()), 
                   paste0("n_",d[,1] %>% unlist()),
                   paste0("NAs_",d[,1] %>% unlist()))

res_num_df <- res %>% as.data.frame()

xlsx::write.xlsx(res_num_df, 
                 "results/Clinical_anlysis/CoMM_dataset/tables/clinical_numeric_variables_in_group3_CoMMpass.xlsx",
                 row.names = F)




#============ categorical ============


library(fastDummies)

data_var_cat <- data %>% select(where(function(x) !is.na(x) && length(unique(x))<15 && length(unique(x))>1), 
                                !!as.symbol(groupvar), 
                                !!as.symbol(ID))

data_var_cat %>% names
data_var_cat_def <- data_var_cat %>% select(# !starts_with("D_TRI") & 
  !ends_with("DETECTED") &
    # !starts_with("HD") & 
    # !starts_with("DEL") & 
    !starts_with("D_IM"),
  # !starts_with("D_QOL"),
  !!groupvar,
) %>% select(- PUBLIC_ID)

data_var_cat_def_dummy <- data_var_cat_def %>% dummy_cols(ignore_na = T)

i <- "SEX_F"

res <- NULL

for(i in data_var_cat_def_dummy %>% names) {
  
  print(i)
  
  try(pval <- with(data_var_cat_def_dummy, {fisher.test(get(groupvar), get(i) )})$p.value %>% round(4))
  
  try(GM <- with(data_var_cat_def_dummy, {gmodels::CrossTable(get(groupvar), get(i)) }) ) 
  
  try(df <- GM$t %>% as.data.frame.matrix() %>% t %>% as.data.frame() %>% rownames_to_column(var = "type")) 
  
  try(freq <- GM$prop.row %>% as.data.frame.matrix() %>% t %>% as.data.frame() %>% set_names(paste0(names(.), "_freq")) ) 
  
  try(tot <- GM$t %>% as.data.frame.matrix() %>% t %>% as.data.frame() %>% rowSums()) 
  
  try(all <- tot %>% sum) 
  
  try(tot_freq <- round(tot/all, 4)) 
  
  try(name_df <- cbind(variable=i, df, p.value=pval, freq, tot, tot_freq))
  
  try(res <- rbind(res, name_df) )
  
}

res_cat_df <- res %>% as.data.frame()

xlsx::write.xlsx(res_cat_df, 
                 "results/Clinical_anlysis/CoMM_dataset/tables/clinical_categoric_variables_in_group3_CoMMpass.xlsx",
                 row.names = F)
