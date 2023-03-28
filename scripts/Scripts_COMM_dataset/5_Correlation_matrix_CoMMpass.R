##################### CoMMpass Correlation matrix SCRIPT ###########################

library(data.table)
library(tidyverse)

# import<-fread("D:/analisi_in_corso/risk_1q_13/CoMMpass/complete_database_1q_13_CoMMpass.txt")
import<-fread("results/complete_database_CoMMpass_1q_13.txt")


# select only relevant variables to associate
names(import)
data<-import %>% select(contains("maj_broad"), contains("maj_focal"), contains("maj_call"), HyperDiploidy, matches("^SeqWGS_"))
names(data)

#____ Build CALLS for traslocations ______

data$t_others <- data$SeqWGS_CCND3_CALL + data$SeqWGS_MAF_CALL + data$SeqWGS_MAFB_CALL

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
                            SeqWGS_CCND3_CALL,
                            SeqWGS_MAF_CALL,
                            SeqWGS_MAFB_CALL
                            ))
names(data)

# relevant calls
data <- as.data.frame(data)
data %>% summarise_all(sum, na.rm =T)
data %>% summarise_all(sum, na.rm =T) %>% as.numeric %>% sort %>% plot + abline(h = 25) + abline(h = 10, col="red")
data %>% summarise_all(sum, na.rm =T) %>% `>`(25) %>% as.vector -> rel
table(rel)

data <- data[rel]
names(data)
data %>% summarise_all(sum, na.rm =T) %>% `>`(25) %>% as.vector %>% table


# add IgH translocated group
data$t_IgH <- data$SeqWGS_WHSC1_CALL + data$SeqWGS_CCND1_CALL + data$t_others

# ------------------ OPTIONAL: add group risk 1q & 13 ------------------------
MM_risk <- ifelse( data$AMP_maj_call_chr_1q == 1 & data$DEL_maj_call_chr_13q ==1, 1,
                   ifelse(data$AMP_maj_call_chr_1q == 1 | data$DEL_maj_call_chr_13q ==1, 2, 3) )

table(MM_risk)
table(MM_risk, data$SeqWGS_WHSC1_CALL)

data$MM_risk_1 <- ifelse(MM_risk==1, 1,0)
data$MM_risk_2 <- ifelse(MM_risk==2, 1,0)
data$MM_risk_3 <- ifelse(MM_risk==3, 1,0)

# dataprep COLUMNS

names(data) <- names(data) %>% str_replace("maj_broad_chr_|maj_call_chr_","") 


asd <- data %>% select(AMP_1p, AMP_1q, AMP_2p:AMP_7q, AMP_8q, AMP_9p,AMP_9q, AMP_11p:AMP_21q, 
                       DEL_1p, DEL_4p:DEL_12q, DEL_13q, DEL_14q, DEL_16p, DEL_16q, DEL_17p, DEL_18p:DEL_22q, 
                       SeqWGS_WHSC1_CALL, SeqWGS_CCND1_CALL, t_others, SeqWGS_MYC_CALL,
                       HyperDiploidy, t_IgH, MM_risk_1:MM_risk_3)


asd <- data %>% select(AMP_1p, AMP_1q, AMP_2p:AMP_7q, AMP_8q, AMP_9p,AMP_9q, AMP_11p:AMP_21q, 
                       DEL_1p, DEL_4p:DEL_12q, DEL_13q, DEL_14q, DEL_16p, DEL_16q, DEL_17p, DEL_18p:DEL_22q, 
                       t_WHSC1=SeqWGS_WHSC1_CALL, 
                       t_CCND1=SeqWGS_CCND1_CALL, 
                       t_others, 
                       t_MYC=SeqWGS_MYC_CALL,
                       HyperDiploidy,
                       t_IgH)


############################# CORRELATION MATRIX ##########################

dataM <- as.matrix(asd)

rownames(dataM) <- rownames(import)


sumsEvents <- apply(dataM, 2, sum, na.rm = TRUE)

0.05*840
42/840

idxEvents <- which(sumsEvents >= 42)

dataMF <- dataM[,idxEvents]
apply(dataMF, 2, sum, na.rm = TRUE)



library(corrplot)

M <- cor(dataMF, use = "pairwise.complete.obs")

pvalueMat <- cor.mtest(M, conf.level = 0.95)

# DEF

pdf("plots/correlation_plots/cor_matrix_CoMM.pdf", 
    width = 10, height = 11)

corrplot(M, type="lower", method="ellipse", 
         tl.col="black", 
         col = colbwr(200),
         p.mat = pvalueMat$p, 
         sig.level = 0.05,
         insig = "pch", 
         pch.cex	=0.8, 
         pch.col ="grey40",
         diag = F,
         tl.srt= 60,
         tl.cex=0.7
         # tl.col=c(rep("orangered",27),rep("deepskyblue2",15),rep("darkorchid3",4), rep("black",5))
         # addgrid.col="white"
) 

dev.off()

