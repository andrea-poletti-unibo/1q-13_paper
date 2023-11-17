#####################  Correlation matrix SCRIPT ###########################

library(data.table)
library(tidyverse)
library(corrplot)


db<-fread("results/complete_database_BO_1q_13.txt")

db_var <- db %>% select(SNP_FILE_NAME, HyperDiploidy, t_4_14, t_11_14, t_14_16, t_14_20, t_6_14)


#____ Build CALLS for traslocations ______
db_var$t_others <- db_var$t_14_16 + db_var$t_14_20 + db_var$t_6_14

db_var$t_IgH <- ifelse(db_var$t_11_14 == 1 | db_var$t_4_14 == 1 | db_var$t_others == 1, 1,0 )


#____ Build CALLS for amp and dels _______

import<-fread("workfiles/BROAD_CALLS/Broad_clonality_calls_BO.txt")

# modify all values in import, assign 1 if > 1 and 0 if < 0.5, if between 0.5 and 1 dont modify
amps <- apply(import %>% select(-sample), c(1,2), function(x) ifelse(x > 1, 1, ifelse(x < 0, 0, x))) %>% as.data.frame()
dels <- apply(import %>% select(-sample), c(1,2), function(x) ifelse(x < -1, 1, ifelse(x > 0, 0, abs(x)))) %>% as.data.frame()

names(amps) <- paste0("Amp_", names(amps))
names(dels) <- paste0("Del_", names(dels))

calls <- cbind(amps, dels)

calls$sample <- import$sample


#_____ add HD and t-IgH to calls ______
db_var <- db_var %>% mutate(sample = SNP_FILE_NAME %>% str_replace_all("\\(|\\)","\\.") %>% str_replace("^([0-9])","X\\1"))
(db_var$sample %in% import$sample) %>% table


allcalls <- left_join(db_var, calls, by = c("sample" = "sample")) %>% select(-SNP_FILE_NAME)





#____ KEEP ONLY UNIQUE/RELVANT CALLS____ 
# unique calls
data <- allcalls %>% select(-c(
    sample,
    t_14_16,
    t_6_14,
    t_14_20
))

names(data)

# relevant calls
data %>% summarise_all(.funs=function(x) sum(x>0.1, na.rm =T)) %>% as.numeric %>% sort %>% plot + abline(h = 25) + abline(h = 10, col="red")
data %>% summarise_all(.funs=function(x) sum(x>0.1, na.rm =T)) %>% ">"(25) %>% as.vector -> rel

25/512
rel %>% table

names(data)[rel]

data <- as.data.frame(data)
data2 <- data[rel]

names(data2)

data2 %>% summarise_all(.funs=function(x) sum(x>0.1, na.rm =T)) %>% `>`(25) %>% as.vector %>% table


#_______ RESHAPE COLUMNS _________
asd <- data2

# reorder columns names
asd <- asd %>% select(Amp_chr_1q:Del_chr_22q,t_4_14,t_11_14,t_others,t_IgH,HyperDiploidy)

# rename columns
names(asd) <- names(asd) %>% str_replace_all("_"," ") 
names(asd) <- names(asd) %>% str_replace_all("chr ","") 

names(asd)




###################################### CORRELATION MATRIX ############################################

dataM <- as.matrix(asd)

rownames(dataM) <- allcalls$sample
rownames(dataM)

sumsEvents <- apply(dataM, 2, function(x) sum(x>0.1, na.rm = TRUE))

table(sumsEvents > 25)
0.05*513
25/513

idxEvents <- which(sumsEvents >= 25)

dataMF <- dataM[,idxEvents]

apply(dataMF, 2, function(x) sum(x>0.1, na.rm = TRUE))



M <- cor(dataMF, use = "pairwise.complete.obs")

pvalueMat <- cor.mtest(M, conf.level = 0.95)
pvalueMat$p

colbwr <- colorRampPalette(c("blue", "white", "red"))
colbwr2 <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                              "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                              "#4393C3", "#2166AC", "#053061") %>% rev) 



#============== create correlation matrix plot ===============

dir.create("plots/correlation_plots/")

M_df <- M %>% as.data.frame(row.names = row.names(M)) %>% rownames_to_column("Alteration")

write_tsv(M_df, "plots/correlation_plots/cor_matrix_clonality_BO.txt")

pdf("plots/correlation_plots/cor_matrix_clonality_BO.pdf", 
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



svg("plots/correlation_plots/cor_matrix_clonality_BO.svg", 
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


 
