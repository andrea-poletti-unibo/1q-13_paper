#####################  Correlation matrix SCRIPT ###########################

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
# MM_risk <- ifelse( data$AMP_maj_call_chr_1q == 1 & data$DEL_maj_call_chr_13q ==1, 1,
#                    ifelse(data$AMP_maj_call_chr_1q == 1 | data$DEL_maj_call_chr_13q ==1, 2, 3) ) 
# 
# table(MM_risk)
# table(MM_risk, data$t_4_14)
# 
# data$MM_risk_1 <- ifelse(MM_risk==1, 1,0)
# data$MM_risk_2 <- ifelse(MM_risk==2, 1,0)
# data$MM_risk_3 <- ifelse(MM_risk==3, 1,0)

# RESHAPE COLUMNS

names(data) <- names(data) %>% str_replace("maj_broad_chr_|maj_call_chr_","") 

asd <- data %>% select(AMP_1q, AMP_3p:AMP_7q, AMP_8q, AMP_9p,AMP_9q, AMP_11p:AMP_21q, 
                       DEL_1p, DEL_6q:DEL_12p, DEL_13q, DEL_14q, DEL_16p, DEL_16q, DEL_17p, DEL_18p, DEL_20p, DEL_22q, 
                       t_4_14, t_11_14, t_others, 
                       HyperDiploidy, t_IgH, MM_risk_1:MM_risk_3)


asd <- data %>% select(AMP_1q, AMP_3p:AMP_7q, AMP_8q, AMP_9p,AMP_9q, AMP_11p:AMP_21q, 
                       DEL_1p, DEL_6q:DEL_12p, DEL_13q, DEL_14q, DEL_16p, DEL_16q, DEL_17p, DEL_18p, DEL_20p, DEL_22q, 
                       t_4_14, t_11_14, t_others, 
                       HyperDiploidy, t_IgH)



############################################### CORRELATION MATRIX #########################################################


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



library(corrplot)

M <- cor(dataMF, use = "pairwise.complete.obs")

pvalueMat <- cor.mtest(M, conf.level = 0.95)
pvalueMat$p

colbwr <- colorRampPalette(c("blue", "white", "red"))
colbwr2 <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                              "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                              "#4393C3", "#2166AC", "#053061") %>% rev) 



# #__________________ 0 ______________
# corrplot(M, type="lower", method="color", tl.col="black", col=colbwr(200), addgrid.col="black") 
# #___________________________________
# 
# 
# corrplot(M, type="lower", addCoef.col = "grey", tl.col="black", number.cex=0.3) # plot correlation matrix


#___________________ 1 _____________

# DEF

pdf("C:/Users/mm_gr/Alma Mater Studiorum Università di Bologna/PROJECT 1q & 13 - Documenti/PAPER_FIGURES_OFFICIAL/4.1b_correlation_matrix_bologna.pdf", 
    width = 10, height = 11)

corrplot(M, type="lower", 
         method="ellipse", 
         tl.col="black", 
         col = colbwr(200),
         p.mat = pvalueMat$p, 
         sig.level = 0.05,
         insig = "pch", 
         pch.cex	=0.8, 
         pch.col ="grey40",
         diag = FALSE,
         tl.srt= 60,
         tl.cex=0.7,
         tl.pos="ld"
         # addgrid.col="white"
         ) 

dev.off()


# ___________________________________
# 
# 
# corrplot(M, type="lower", method="number", tl.col="black") 
# corrplot(M, type="lower", method="shade", tl.col="black") # plot correlation matrix
# 
# 
# #___________________ 2 _____________
# corrplot(M, type="lower", method="color", tl.col="black", # addCoef.col = "grey40", number.cex=0.45, 
#          col=colbwr2(200),
#          p.mat = pvalueMat$p, 
#          sig.level = 0.05,
#          insig = "pch", pch.cex	=1.2,
#          diag = FALSE,
#          addgrid.col="black") 
# #___________________________________
# 
# 
# corrplot(M, type="lower", method="pie", tl.col="black") 
# 
# corrplot(M, type="lower", method="pie", outline=T, addgrid.col="white" , tl.col="black") 
# 
# corrplot(M, type="lower", order="hclust", tl.col="black") 
# corrplot(M, order="hclust", addrect = 6, tl.col="black") 
# 
# 
# #_________________ 3 __________________
# corrplot(M, type="lower", method="pie", outline=T, 
#          col=colbwr2(200),
#          p.mat = pvalueMat$p, 
#          sig.level = 0.05,
#          insig = "blank",
#          tl.col="black",
#          diag = FALSE
# ) 
# #______________________________________
# 
# 
# 
# 
# 
# corrplot(M2, type="lower") # plot CLUSTERD correlation matrix
# 
# # other plots (mixed)
# corrplot.mixed(M2)
# corrplot.mixed(M2, lower="ellipse", upper="circle")
# corrplot.mixed(M2, lower="square", upper="circle")
# corrplot.mixed(M2, lower="shade", upper="circle")
# corrplot.mixed(M2, tl.pos="lt")
# corrplot.mixed(M2, tl.pos="lt", diag="u")
# corrplot.mixed(M2, tl.pos="lt", diag="l")
# corrplot.mixed(M2, tl.pos="n")
# 
# 
