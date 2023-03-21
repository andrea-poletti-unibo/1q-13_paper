#####################  FISHER matrix SCRIPT ###########################

library(data.table)
library(tidyverse)

# import<-fread("D:/analisi_in_corso/risk_1q_13/complete_database_1q_13_181019.txt")
import<-fread("C:/Users/mm_gr/Alma Mater Studiorum Università di Bologna/PROJECT 1q & 13 - Documenti/complete_database_1q_13_181019.txt")


# remove spaces from SNP names 
import$nSNP_1q_13_project
rownames(import) <- import$nSNP_1q_13_project %>% str_replace(" ", "_")
rownames(import)

# select only relevant variables to associate
names(import)
data<-import %>% select(contains("maj_broad"), contains("maj_focal"), HyperDiploidy, matches("^t_"))
names(data)

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

#____ Build CALLS for traslocations ______

data$t_others <- data$t_14_16 + data$t_14_20 + data$t_6_14

#____ KEEP ONLY UNIQUE/RELVANT CALLS____ 

# unique calls
data <- data %>% select( -c(AMP_maj_broad_chr_1q, 
                            AMP_maj_broad_chr_8q, 
                            DEL_maj_broad_chr_17p, 
                            DEL_maj_broad_chr_13q, 
                            DEL_maj_broad_chr_14q, 
                            DEL_maj_broad_chr_1p,
                            DEL_maj_broad_chr_16q, 
                            contains("maj_focal"),
                            t_14_16,
                            t_6_14,
                            t_14_20
                            ))
names(data)

# relevant calls
data <- as.data.frame(data)
data %>% summarise_all(sum, na.rm =T) %>% as.numeric %>% sort %>% plot + abline(h = 25) + abline(h = 10, col="red")
data %>% summarise_all(sum, na.rm =T) %>% `>`(25) %>% as.vector -> rel

data <- data[rel]
names(data)
data %>% summarise_all(sum, na.rm =T) %>% `>`(25) %>% as.vector %>% table


#------------------ OPTIONAL: add group risk 1q & 13 ------------------------
MM_risk <- ifelse( data$AMP_maj_call_chr_1q == 1 & data$DEL_maj_call_chr_13q ==1, 1,
                        ifelse(data$AMP_maj_call_chr_1q == 1 | data$DEL_maj_call_chr_13q ==1, 2, 3) ) 

table(MM_risk)
table(MM_risk, data$t_4_14)

data$MM_risk_1 <- ifelse(MM_risk==1, 1,0)
data$MM_risk_2 <- ifelse(MM_risk==2, 1,0)
data$MM_risk_3 <- ifelse(MM_risk==3, 1,0)

# RESHAPE COLUMNS

names(data) <- names(data) %>% str_replace("maj_broad_chr_|maj_call_chr_","") 

#______ EXPLORATIVE > all variables ______
asd <- data %>% select(AMP_1q, AMP_3p:AMP_7q, AMP_8q, AMP_9p,AMP_9q, AMP_11p:AMP_21q, 
                       DEL_1p, DEL_6q:DEL_12p, DEL_13q, DEL_14q, DEL_16p, DEL_16q, DEL_17p, DEL_18p, DEL_20p, DEL_22q, 
                       t_4_14, t_11_14, t_others, 
                       HyperDiploidy, t_IgH, MM_risk_1:MM_risk_3)

#______ PAPER > no MMrisk groups ______

asd <- data %>% select(AMP_1q, AMP_3p:AMP_7q, AMP_8q, AMP_9p,AMP_9q, AMP_11p:AMP_21q, 
                       DEL_1p, DEL_6q:DEL_12p, DEL_13q, DEL_14q, DEL_16p, DEL_16q, DEL_17p, DEL_18p, DEL_20p, DEL_22q, 
                       t_4_14, t_11_14, t_others, 
                       HyperDiploidy, t_IgH)



# Fisher MATRIX
dataM <- as.matrix(asd)
rownames(dataM) <- rownames(import)
rownames(dataM)

sumsEvents <- apply(dataM, 2, sum, na.rm = TRUE)
sumsEvents

0.05*513
25/513

idxEvents <- which(sumsEvents >= 25)

dataMF <- dataM[,idxEvents]
apply(dataMF, 2, sum, na.rm = TRUE)

# table(dataMF[,1], dataMF[,3])
# f <- fisher.test(table(dataMF[,1], dataMF[,2]))
# c <- chisq.test(table(dataMF[,1], dataMF[,2]))
# 
# f$p.value
# c$p.value


i=1
j=3

mat <- matrix(data = NA, nrow= ncol(dataMF), ncol = ncol(dataMF), dimnames = list(colnames(dataMF),colnames(dataMF)))

for (i in 1:ncol(dataMF)){
  for(j in 1:ncol(dataMF)){
    f <- fisher.test(dataMF[,i], dataMF[,j], alternative = "greater")
    p <- f$p.value
    mat[i,j] <- p
  }
}


library(reshape2)

mat.UpTri<- mat
mat.UpTri[upper.tri(mat.UpTri, diag = T)] <- NA

mat.melt <- as.data.frame(melt(mat.UpTri))

mat.melt$logP <- log(mat.melt$value)

# create a thresholded version of p-value to handle the extremely significative values
mat.melt$pval.thresh <- mat.melt$value
mat.melt$pval.thresh[mat.melt$pval.thresh<0.001] <- 0.001

mat.melt$Var_fact <- mat.melt$Var1 %>% factor(levels=(mat.melt$Var1 %>% rev %>% unique))

ggplot(mat.melt, aes(mat.melt$Var_fact, mat.melt$Var2, fill=mat.melt$pval.thresh))+
  geom_tile(color="white")+
  scale_fill_gradient2(low="red", high="black", mid="orange", limits= c(0,0.05), midpoint = 0.025, space="Lab", name="Fisher p-value", na.value = "white")+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90, vjust=0,hjust=1))+
  coord_fixed()+
  labs(x="", y="") 



#=============  fisher exclusive assocition (inverted test) =================

mat_inv <- matrix(data = NA, nrow= ncol(dataMF), ncol = ncol(dataMF), dimnames = list(colnames(dataMF),colnames(dataMF)))
for (i in 1:ncol(dataMF)){
  for(j in 1:ncol(dataMF)){
    f <- fisher.test(dataMF[,i], dataMF[,j], alternative = "less")
    p <- f$p.value
    mat_inv[i,j] <- p
  }
}

mat.UpTri_inv<- mat_inv
mat.UpTri_inv[upper.tri(mat.UpTri_inv, diag = T)] <- NA

mat.melt_inv <- as.data.frame(melt(mat.UpTri_inv))

mat.melt_inv$logP <- log(mat.melt_inv$value)

mat.melt_inv$pval.thresh <- mat.melt_inv$value
mat.melt_inv$pval.thresh[mat.melt_inv$pval.thresh<0.001] <- 0.001


mat.melt_inv$Var_fact <- mat.melt_inv$Var1 %>% factor(levels=(mat.melt_inv$Var1 %>% rev %>% unique))

ggplot(mat.melt_inv, aes(mat.melt_inv$Var_fact, mat.melt_inv$Var2, fill=mat.melt_inv$pval.thresh))+
  geom_tile(color="white")+
  scale_fill_gradient2(low="blue", high="purple", mid="skyblue2", limits= c(0,0.05), midpoint = 0.025, space="Lab", name="Fisher p-value", na.value = "white")+
  theme_bw()+
  theme(axis.text.x=element_text(angle=60, vjust=1, size=10, hjust=1))+
  coord_fixed()+
  labs(x="", y="")



source("D:/Cartelle_personali/Andrea/GIT_Shirke019/Various_analysis/new_aes_multiple_colorscales_ggplot.R")
# source("C:/Users/andre/Desktop/MM/R/PROJECTS/MM_workspace/Various_analysis/new_aes_multiple_colorscales_ggplot.R")


#===== JOINT PLOT ========

mat.melt$pval.thresh_inv <- mat.melt_inv$pval.thresh

ggplot(mat.melt, aes(Var_fact, Var2))+
  theme_bw()+
  
  geom_tile(color="white", aes(fill=pval.thresh), data= mat.melt[mat.melt$pval.thresh<0.05 | is.na(mat.melt$pval.thresh_inv),])+
  scale_fill_gradient2(low="red", high="black", mid="orange", limits= c(0,0.05), midpoint = 0.025, space="Lab", name="Fisher p-value", na.value = "white")+

  new_scale("fill") +
  
  geom_tile(color="white", aes(fill= pval.thresh_inv),  data = mat.melt[mat.melt$pval.thresh_inv <0.05 | is.na(mat.melt$pval.thresh_inv) ,])+
  scale_fill_gradient2(low="blue", high="purple", mid="skyblue2", limits= c(0,0.05), midpoint = 0.025, space="Lab", name="Fisher p-value", na.value = "white")+
  
  theme(axis.text.x=element_text(angle=60, vjust=1, size=10, hjust=1),
        panel.border = element_blank(),
        legend.position = c(0.9, 0.8))+
  coord_fixed()+
  labs(x="", y="")


ggsave("C:/Users/mm_gr/Alma Mater Studiorum Università di Bologna/PROJECT 1q & 13 - Documenti/PAPER_FIGURES_OFFICIAL/4.1_fisher_matrix_bologna.pdf", 
       dpi = 300, width = 10, height = 10, units = "in")





#################### Fisher boshloo test ################################
# 
# library(exact2x2)
# 
# 
# t<-table(data$AMP_broad_chr_1q, data$t_4_14)
# 
# t
# t[1,2]
# 
# 
# m <- matrix(c(t[1,1], t[1,1]+t[1,2], t[2,1], t[2,1]+t[2,2]),
#        nrow = 2)
# 
# m
# 
# mt <- matrix(c(t[1,1], t[2,1], t[1,2], t[2,2]),
#              nrow = 2)
# mt
# t
# 
# fisher.test(m, alternative = "greater")
# fisher.test(t, alternative = "greater" )
# 
# fisher.test(data$AMP_broad_chr_1q, data$t_4_14)
# 
# 
# m
# fisher.test(m, alternative = "greater" )
# 
# 
# fisher_test <- fisher.test
# 
# bosch_test <- boschloo(t[1,1], t[1,1]+t[1,2], t[2,1], t[2,1]+t[2,2]) 
# 
# bosch_test
# 
# t
# t[1,2]
# 
# 
# library(Exact)
# 
# exact.test(t, alternative = "greater", method="Boschloo", npNumbers=1000)
# 
# 
# 
# data <- matrix(c(7, 8, 12, 3), 2, 2, byrow=TRUE)
# exact.test(data, alternative="less",to.plot=TRUE)
# exact.test(data, alternative="two.sided", interval=TRUE, beta=0.001, npNumbers=100,
#            method="Z-pooled",to.plot=FALSE)
# exact.test(data, alternative="two.sided", interval=TRUE, beta=0.001, npNumbers=100,
#            method="Boschloo", to.plot=FALSE)
# 
# fisher.test(data, alternative = "two.sided")
