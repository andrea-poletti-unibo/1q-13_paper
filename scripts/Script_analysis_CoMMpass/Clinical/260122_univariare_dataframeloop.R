#######################
# COMMPASS - UNIVARIATE - TABLE
# ANALYSIS OF 260122 
##########################
library(sjPlot)

#install.packages("sjPlot")
library(stargazer)


#set directory of the project
getwd()
#ajsi
setwd("C:/Users/ajsik/Desktop/University (non online)/Tirocinio Sant'Orsola (non online)/git_mmrisk")


#create temp functions for COX univariate
tmpfun_os <- function(x) as.formula(paste("OS",x,sep="~"))
tmpfun_pfs <- function(x) as.formula(paste("PFS",x,sep="~"))


#import the CoMMpass dataset with all variables
load(file="./data/co2701.Rda")
OS_co <- Surv(import$OS_months, import$censos)
PFS_co <- Surv(import$PFS_months, import$censpfs)


#import the Bologna dataset with all variables
load(file="./data/bo2701.Rda")
OS_bo <- Surv(import$OS_MONTHS, import$OS_EVENT)
PFS_bo <- Surv(import$PFS_I_MONTHS, import$PFS_I_EVENT)


#1.___create string with the list of variables_________
OS_co->OS
PFS_co->PFS

#vars of interest of the univariate analysis
var_co<-c("OLD",
       "SEX",
       "D_LAB_serum_beta2_microglobulin", 
       "D_LAB_chem_albumin",
       "HB_m_105",
       "PLT_m_150",
       "PC_M_60",
       "D_LAB_chem_ldh",
       "R_ISS",
       "D_IM_LIGHT_CHAIN_BY_FLOW2",
       "sctflag",
       "NOTE_TRANSLOCATION",
       "Hyperdiploid",
       "`DEL_maj_broad_chr_17p`",
       "`DEL_maj_broad_chr_1p`",
       "`DEL_maj_focal_TP53`",
       "`AMP_maj_broad_chr_1q`",
       "`DEL_maj_broad_chr_13q`",
       "MMrisk_CLASSIFIER_ALL")



#_______txt output___________initialize output____________
with(import, coxph(OS ~ OLD )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .) %>% cbind(N=sum(!is.na(import$OLD))) %>% cbind(variable="OLD") %>% write_tsv(paste0(outpath,"report_cox_analysis_260122.txt"), append = T, col_names = T)
with(import, coxph(PFS ~ OLD )) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .) %>% cbind(N=sum(!is.na(import$OLD))) %>% cbind(variable="OLD") %>% write_tsv(paste0(outpath,"report_cox_analysis_260122.txt"), append = T, col_names = F)

i="OLD"

for (i in var_co) {
  print(i)
  coxph(tmpfun_os(i), data = co)%>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .) %>% cbind(N=sum(!is.na(i))) %>% cbind(variable=i) %>% write_tsv(paste0(outpath,"report_cox_analysis_260122.txt"), append = T, col_names = F)
  coxph(tmpfun_pfs(i), data = co) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .) %>% cbind(N=sum(!is.na(i))) %>% cbind(variable=i) %>% write_tsv(paste0(outpath,"report_cox_analysis_260122.txt"), append = T, col_names = F)
  
}


#________dataframe output_____________
output1<-output2<-NULL
for (i in var_co) {
  print(i)
  out1<-coxph(tmpfun_os(i), data = co) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .) 
  out2<-coxph(tmpfun_pfs(i), data = co) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)
  output1<-rbind(out1, out2)
  output2<-rbind(output1, output2)
}
output2
output<-output2[,c(1,2,3,7,8,6)]
names(output)<-c("OS/PFS", "variable", "HR","Lower CI", "Upper CI", "p-value")
output_co<-output
output_co


save(file="./Script_analysis_CoMMpass/Clinical/Univariate_Rda_CoMMpass/univariate_co.Rda", output_co)

######################

#_________latex output_________
xtable::xtable(output_co)




##_________Bologna___________
OS_bo->OS
PFS_bo->PFS

#1.__create variables____ of interest of univariate analysis
var_bo<-c("PROTOCOL_REV",
        "OLD",
        "SEX",
        "B2M",
        "ALBUMIN",
        "HB_m_105",
        "PLT_m_150",
        "PC_M_60",
        "LDH_UL",
        "ISS",
        "R_ISS",
        "IG_ISOTYPE_REV",
        "LIGHT_CHAIN",
        "ASCT",
        "ASCT_MMrisk",
        "NUMBERS_TX_FRONT_LINE",
        "MAINTENANCE",
        "CONSOLIDATION",
        "Induction_therapy",
        "Induction_response_M_VGPR",
        "NOTE_TRANSLOCATION",
        "`DEL_maj-broad_chr_17p`",
        "`DEL_maj-focal_TP53`",
        "FISH_Del_17p",
        "FISH_Del_1p36",
        "`AMP_maj-broad_chr_1q`",
        "`DEL_maj-broad_chr_13q`",
        "MMrisk_class",
        "MMrisk_AB_ALL")
names(bo)        
i="PROTOCOL_REV"

output1<-output2<-NULL
for (i in var_bo) {
  print(i)
  out1<-coxph(tmpfun_os(i), data = bo) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .) 
  out2<-coxph(tmpfun_pfs(i), data = bo) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)
  output1<-rbind(out1, out2)
  output2<-rbind(output1, output2)
}
output2
output<-output2[,c(1,2,3,7,8,6)]
names(output)<-c("OS/PFS", "variable", "HR","Lower CI", "Upper CI", "p-value")
output_bo<-output
output_bo$`p-value`<-round(output$`p-value`,3)

#save(file="./Script/Clinical/Univariate_Rda_Bologna/univariate_bo.Rda", output_bo)

xtable::xtable(output_bo,digits = 3)
sjPlot::sjt.df(output_bo)
sjPlot::tab_df(output_bo)

levels(bo$MMrisk_AB_ALL)
