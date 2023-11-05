#####################  Correlation matrix SCRIPT ###########################

library(data.table)
library(tidyverse)

# import<-fread("D:/analisi_in_corso/risk_1q_13/complete_database_1q_13_181019.txt")
import<-fread("results/complete_database_BO_1q_13.txt")


# select only relevant variables to associate
names(import)
data<-import %>% select(contains("maj_broad"), contains("maj_focal"), contains("maj_call"), HyperDiploidy, matches("^t_"))
names(data)

#____ Build CALLS for traslocations ______

data$t_others <- data$t_14_16 + data$t_14_20 + data$t_6_14

data$t_IgH <- ifelse(data$t_11_14 == 1 | data$t_4_14 == 1 | data$t_others == 1, 1,0 )

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

asd <- data %>% select(AMP_1q, AMP_3p:AMP_7q, AMP_8q, AMP_9p,AMP_9q, AMP_11p:AMP_21q, 
                       DEL_1p, DEL_6q:DEL_12p, DEL_13q, DEL_14q, DEL_16p, DEL_16q, DEL_17p, DEL_18p, DEL_20p, DEL_22q, 
                       t_4_14, t_11_14, t_others, 
                       HyperDiploidy, t_IgH, MM_risk_1:MM_risk_3)


asd <- data %>% select(AMP_1q, AMP_3p:AMP_7q, AMP_8q, AMP_9p,AMP_9q, AMP_11p:AMP_21q, 
                       DEL_1p, DEL_6q:DEL_12p, DEL_13q, DEL_14q, DEL_16p, DEL_16q, DEL_17p, DEL_18p, DEL_20p, DEL_22q, 
                       t_4_14, t_11_14, t_others, 
                       HyperDiploidy, t_IgH)


names(asd) <- names(asd) %>% str_replace("AMP","Amp") 
names(asd) <- names(asd) %>% str_replace("DEL","Del") 
names(asd) <- names(asd) %>% str_replace("_"," ") 
names(asd) <- names(asd) %>% str_replace("_",";") 
names(asd)

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



#============== create correlation matrix plot ===============

dir.create("plots/correlation_plots/")

pdf("plots/correlation_plots/cor_matrix_BO.pdf", 
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



svg("plots/correlation_plots/cor_matrix_BO.svg", 
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



