
# load packages
library(data.table)
library(tidyverse)
library(survival)
library(survminer)
library(data.table)
library(readr)
library(edgeR)
library(limma)
library(Glimma)
library(RColorBrewer)
library(biomaRt)
library(knitr)
library(lattice)


# load database
import <- fread("results/complete_database_CoMMpass_1q_13.txt")

# define the MM risk label from commpass broad + focal CNAs
import$AMP_1q_genes_all <- with(import, 
                                ifelse(`AMP_maj_focal_ANP32E`==1 | `AMP_maj_focal_MCL1` ==1 | `AMP_maj_focal_CKS1B`==1, 1,0 ))

import$MMrisk_1q_all <- ifelse( import$`AMP_maj_broad_chr_1q` == 1 | import$AMP_1q_genes_all ==1, 1, 0)
import$MMrisk_13_all <- ifelse( import$`DEL_maj_broad_chr_13q` ==1 | import$`DEL_maj_focal_RB1` ==1, 1, 0)

import$MMrisk_CLASSIFIER_ALL <- 3  - import$MMrisk_1q_all - import$MMrisk_13_all
import$MMrisk_CLASSIFIER_ALL <- dplyr::recode(as.character(import$MMrisk_CLASSIFIER_ALL), "1"="risk1","2"="risk2","3"="risk3")



# define the dataframe with group assignments    
compassGroup <- import %>% dplyr::select(Study_Visit_iD,MMrisk_CLASSIFIER_ALL) 

NAME<-"CoMMpass_risk1_vs_risk3" # set name of analysis

path <- "results/DEG_analysis/"

dir.create(path)

opts_knit$set(root.dir = path)

DIR<- path

setwd(path)



# load GENE COUNT trascriptome data from CoMMpass RNA-seq
importCounts <- fread("../../input_data/CoMMpass/MMRF_CoMMpass_IA13a_E74GTF_Salmon_Gene_Counts.txt.gz") # import gene counts data

# define the dataframe with group assignments
groupImport<- compassGroup

# 1.2 ____REORDERING SAMPLE COLUMNS - CRITICAL STEP_____

importCounts<-as.data.frame(importCounts)
head(names(importCounts))

ord<-order(names(importCounts))

importCounts<- importCounts[,ord] # reordering step!

head(names(importCounts))





# 2.1  integration of ENSEMBL gene annotation with HGNC_SYMBOL, START, END info downloaded from bioMart 

genelist<-importCounts$GENE_ID #create commpass gene list

human <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

gene_coords=getBM(attributes=c("hgnc_symbol","ensembl_gene_id", "start_position","end_position"),
                  filters="ensembl_gene_id",
                  values=genelist, #download the additional annotation for the compass gene list
                  mart=human)

#import "bioMart_genes_annot.RData" prepared file if no internet connection available 
# load(paste0(DIR,"/bioMart_genes_annot.RData")) #NOTA BENE!!!!!

#creation of length variable
gene_coords$length = gene_coords$end_position - gene_coords$start_position

#rename the gene header of the import file
names(importCounts)[1]<-"ensembl_gene_id"

#merge the import file with the new downloaded annotations
mergedInfoCounts<- left_join(gene_coords, importCounts, by = "ensembl_gene_id") # NB! FROM 57997 genes TO 52460 genes





# 2.2  Select only DIAGNOSIS samples from CoMMpass 
onlyDiagnosis<- dplyr::select(mergedInfoCounts, contains("_1_BM")) 

samplenamesD<-names(onlyDiagnosis)

prepData<- onlyDiagnosis # CREATION of prepData object




# 2.2b _______ [OPTIONAL] Exclude samples with translocations involving cycline deregulations ________
# t_import <- import %>% filter(SeqWGS_WHSC1_CALL == 1 | SeqWGS_CCND1_CALL==1 | SeqWGS_CCND3_CALL ==1 | SeqWGS_MAF_CALL==1 | SeqWGS_MAFB_CALL==1)
# t_import$SeqWGS_WHSC1_CALL + t_import$SeqWGS_CCND1_CALL + t_import$SeqWGS_CCND3_CALL + t_import$SeqWGS_MAF_CALL + t_import$SeqWGS_MAFB_CALL
# 
# drops <- t_import$Study_Visit_iD
# 
# filtered<-prepData[, !(names(prepData) %in% drops)] #exclude translocated patients from prepData
# 
# #checks
# sum(names(filtered)%in% drops)
# sum(names(filtered)%in%names(prepData)==F)
# samplenamesDnoT<-names(filtered)
# 
# prepData<-filtered # !!! SAVE prepData !!!

# # 2.2c _________ [OPTIONAL] others sample type filtering steps ____________
# #define the condition for keeping samples
# keep<- dplyr::filter(groupImport, groupImport$`1q` == 0) %>% dplyr::select(sample) %>% unlist()
# 
# prepData<- prepData[, names(prepData) %in% keep]




# 2.3 ________ define GROUP variable in order to perform differential expression analysis between groups ________

groupSelect<- dplyr::filter(groupImport, Study_Visit_iD %in% names(prepData)) # select the name list to filter


samplewithgroup<- groupSelect$Study_Visit_iD
## groupCheck<-filter(groupImport, sample %in% samplewithgroup)

#excluding samples without a group
onlyGroup<- dplyr::select(prepData, samplewithgroup)

prepData<-cbind(mergedInfoCounts[,c(1,2,5)],onlyGroup) # !!! SAVE prepData !!! + add 3 gene annotations cols 

dim(prepData)# row= genes, cols= info + samples


# 2.4________  create essential variable necessary for DEG analysis ________

#_______essential 1: group
group<-as.factor(groupSelect$MMrisk_CLASSIFIER_ALL) 

length(group)
table(group)

#_______essential 2: counts
counts<-as.matrix(prepData[,c(4:ncol(prepData))])
dim(counts)

#_______essential 3: samples
samples<-data.frame(samples=names(prepData[-c(1:3)]))

#_______essential 4: genes
genes<-data.frame(genes=prepData[,1:3])

# analyze if there are duplicated genes
notdup<-genes[!duplicated(genes$genes.ensembl_gene_id),]
dup<-genes[duplicated(genes$genes.ensembl_gene_id),]



# 3.1 ________ data import in edgeR object________

x<- DGEList(counts,samples=samples,genes=genes)
dim(x)

# group assignment
x$samples$group<- group



# 3.2 ________ data pre-processing ________

# trasformations from the raw scale to Counts Per Million (CPM) and logCMP
cpm<-cpm(x)
lcpm<-cpm(x, log=T)

#notice differences after trasnformation
cpm[1:5,1:5]
lcpm[1:5,1:5]
counts[1:5,1:5]


# #~~~~~~~~~ OPTIONAL: plot a GRID to best choose the optimum parameters to filter lowexpr genes ~~~~~~~~~
# #
# i_cpm<-vector()
# j_pat<-vector()
# c_keep<-vector()
# for (i in seq(0.1,2,by=0.05)){   # try different CPM thresholds
#   for (j in seq(1,30,by=1)){     # try different min patients cutoffs
#     vec <- rowSums(cpm>i)>=j
#     c<-table(vec)[2]             # use the kept genes as the result value
#     c_keep<-append(c_keep,c)
#     i_cpm<-append(i_cpm,i)
#     j_pat<-append(j_pat,j)
#   }
# }
# 
# plot(c_keep)
# 
# #creation of a dataframe
# data <- data.frame(i_cpm, j_pat, c_keep)
# 
# #plot with ggplot
# ggplot(data, aes(i_cpm, j_pat, z= c_keep)) +
#   geom_tile(aes(fill = c_keep)) +
#   theme_bw()
# 
# #plot with lattice
# library(lattice)
# par(mar=c(3,4,2,2))
# levelplot(c_keep ~ i_cpm*j_pat, data=data  , xlab="cpm" , col.regions = heat.colors(100)[length(heat.colors(100)):1]   , main="")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



# 3.3 ________removing low expressed genes________

# see how many genes have 0 counts in every sample 
table(rowSums(x$counts==0)== ncol(x)) # ncol(x) = total number of samples 

# filter genes basing on low-expr thresholds
keep.exprs<- rowSums(cpm>1)>=5 # select the genes that have a cpm >1 (typical CPM threshold) in at least 5 samples (minimum significative group size)
table(keep.exprs)

x<-x[keep.exprs,, keep.lib.sizes=FALSE]
dim(x)

# # plotting the results of filtering
# library(RColorBrewer)
# nsamples <- ncol(x)
# col <- brewer.pal(nsamples, "Paired")
# par(mfrow=c(1,2))
# plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.21), las=2,
#      main="", xlab="")
# title(main="A. Raw data", xlab="Log-cpm")
# abline(v=0, lty=3)
# for (i in 2:nsamples){
#   den <- density(lcpm[,i])
#   lines(den$x, den$y, col=col[i], lwd=2)
# }
# lcpm <- cpm(x, log=TRUE)
# plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.21), las=2,
#      main="", xlab="")
# title(main="B. Filtered data", xlab="Log-cpm")
# abline(v=0, lty=3)
# for (i in 2:nsamples){
#   den <- density(lcpm[,i])
#   lines(den$x, den$y, col=col[i], lwd=2)
# }




# 3.4 ________ Normalization across samples in the experiment ________

x<- calcNormFactors(x, method = "TMM")
head(x$samples$norm.factors, n=20)


# # **PLOT> intersamples expr normalization boxplots**
# # explore the effect of normalization (duplicated and modified data to enhance the visual effect)
# x2<-x[,1:50] # ONLY first 20 samples in this example
# x2$samples$norm.factors <- 1
# # x2$counts[,1]<- ceiling(x2$counts[,1]*0.05) # first sample counts are reduced to 5% of original values
# # x2$counts[,2]<-x2$counts[,2]*5 # second sample counts are inflated to 500% of original value
# 
# nsamples <- ncol(x)
# col <- brewer.pal(nsamples, "Paired")
# 
# par(mfrow=c(1,2))
# lcpm <- cpm(x2, log=TRUE)
# boxplot(lcpm, las=2, col=col, main="")
# title(main="A. Example: Unnormalised data",ylab="Log-cpm")
# 
# x2 <- calcNormFactors(x2)
# x2$samples$norm.factors
# lcpm <- cpm(x2, log=TRUE)
# boxplot(lcpm, las=2, col=col, main="")
# title(main="B. Example: Normalised data",ylab="Log-cpm")


# [4] =================== limma & Glimma: identify DEG and visualize results ==========================

# 4.1________ Multi Dimensional Scaling (MDS) analysis plot (~PCA of the groups)________

library(limma)

# (search initial global evidence for a differential expression between groups)
lcpm <- cpm(x, log=TRUE)


# **PLOT> MDS plot with limma**

# # MDS PLOT

pdf("limma_gene_expression_MDS.pdf", width = 10, height = 5)

par(mfrow=c(1,2))
col.group <- group
levels(col.group) <- brewer.pal(nlevels(col.group), "Set1")
col.group <- as.character(col.group)

plotMDS(lcpm, #labels=group,
        pch=20,
        col=col.group,
        dim=c(1,2))# plot dimension 1 vs dimension 2


plotMDS(lcpm, #labels=group,
        pch=20,
        col=col.group,
        dim=c(2,3)) # plot dimension 2 vs dimension 3

title(main="Limma gene expression MDS plot")

dev.off()

# **PLOT> MDS INTERACTIVE plot with Glimma**

# Glimma interactive MDS viz - ALL dimensions
glMDSPlot(lcpm, labels="group", groups=x$samples$group,
          launch=TRUE) #launch = TRUE opens the plot in a browser



# 4.2________ Differential expression experimental settings: DESIGN and CONTRAST matrices ________

#create a DESIGN MATRIX from group

design<- model.matrix(~0+group)
colnames(design)<-gsub("group","",colnames(design))

design[1:10,]

#create a CONTRAST MATRIX with comparison of interest between groups

# SET the groups

contr.matrix<- makeContrasts(
  r1v3 = `risk1`-`risk3`,
  r1v2 = `risk1`-`risk2`,
  r2vr3 = `risk2`-`risk3`,
  levels = colnames(design))

contr.matrix


# 4.3________ VOOM function: removing HETEROSCEDASCITY from count data (and transform count data in logCPM) ________

v<- voom(x, design, plot=TRUE)
v[1:10,1:10]



# 4.4________Fitting linear models and bayesian t-test for comparisons of interest________

vfit <- lmFit(v, design) # calculate the mean expression level for each gene in the different groups

vfit<- contrasts.fit(vfit, contrasts = contr.matrix) # compare the groups

efit<- eBayes(vfit) # calculate statistical bayesian-variance-adjusted T-test to find significant differential expressed genes (pval<0.05)


# **PLOT> voom plot: median-variance heteroscedascity correction**
plotSA(efit)

# **RESULTS> DE genes**
# Examining the number of DE genes: Significance is defined using an adjusted p-value cutoff that is set at 5%
summary(decideTests(efit))

tab= topTable(efit, coef = 1, number=1000, adjust.method = "BH")
# View(tab)

topGenes<- tab[tab[,"adj.P.Val"]< 0.001,]


# **PLOT> volcano plot of DEG (p<0.05)**
dev.off()
pdf("volcano_500_DEG.pdf", width = 9, height = 7)
volcanoplot(efit,coef = 1, highlight = 500, names=efit$genes$genes.hgnc_symbol)
title(main="Volcano plot of top 500 deregulated genes")
dev.off()



# 4.4b ________Strictier definition of significance ________
# "treat" method to calculate p-values from empirical Bayes moderated t-statistics with a minimum log-FC requirement. 

tfit <- treat(vfit, lfc=1) # lfc=Log2FoldChange (lfc=1 is then equal to 2 fold-change)

# **RESULTS> strict DE genes**

#Examining the number of strict DE genes (with minimum fold-change > 2)
dt <- decideTests(tfit)
summary(dt)

strictNum <- length(which(dt[,1]!=0))

tabStrict<- topTreat(tfit, coef=1, sort.by = "p", number=strictNum) # SET the correct number of strict DEG genes !
# View(tabStrict)


# **PLOT> volcano plot of STRICT DEG**
dev.off()
pdf("volcano_500_strict-DEG.pdf", width = 9, height = 7)
volcanoplot(tfit,coef = 1, highlight = strictNum, names= tfit$genes$genes.hgnc_symbol) # SET the correct number of strict DEG genes
title(main="Volcano plot of STRICT (fold-change>2) deregulated genes")
dev.off()


# 4.5 ________ EXTRACT both upregolated and downregolated STRICT significative DEG ________

# **RESULTS> gene LISTS**

# common up and down DEG
de.common<- which(dt[,1]!=0)
length(de.common)

tfit$genes$genes.hgnc_symbol[de.common]

# up-regulated genes
de.up<- which(dt[,1]==1)
length(de.up)

tfit$genes$genes.hgnc_symbol[de.up]

# down-regulated genes
de.down<- which(dt[,1]== -1)
length(de.down)

tfit$genes$genes.hgnc_symbol[de.down]


# Examining individual DE genes from top to bottom
DEGresults<- topTreat(tfit, coef=1, n=Inf)
# View(DEGresults)



# 4.6 ________ EXPORT the strict DE gene list and the analysis results ________
write.fit(tfit, dt, file=paste0("DEG_results_",NAME,".txt")) # analysis result
write_tsv(tabStrict,paste0("DEG_genes_",NAME,".txt")) # strict DE list

# no traslocation 
# write.fit(tfit, dt, file=paste0("DEG_results_noT_",NAME,".txt")) # analysis result
# write_tsv(tabStrict,paste0("DEG_genes_noT_",NAME,".txt")) # strict DE list


# 4.7 ________Useful graphical representations of differential expression results ________
# **PLOT> MD plot of STRICT DEG by limma**
# MD plot with UPREGOLATED and DOWNREGOLATED strict significative genes
plotMD(tfit, column=1, status=dt[,1], main=colnames(tfit)[1], xlim=c(-8,13))

# **PLOT> MD plot of STRICT DEG by Glimma**
# Glimma interactive MD plot 
glMDPlot(tfit, coef=1, status=dt, main=colnames(tfit)[1],
         side.main="genes.hgnc_symbol", counts=x$counts, groups=group, launch=TRUE) # Glimma INTERACTIVE




# 5.0____________________ SAVE ENV RDATA __________________

# dir.create("../../workfiles/DEG/", recursive = T, showWarnings = F)

# save.image("../../workfiles/DEG/RNA-seq_DEG_coMMpass_env.Rdata")

# save for no Traslocation
# save.image("../../workfiles/DEG/RNA-seq_DEG_coMMpass_env_noT.Rdata")


#________________________ LOAD ENV RDATA _________________
# load DEG pipeline .Rdata
# load("../../workfiles/DEG/RNA-seq_DEG_coMMpass_env.Rdata")

# load DEG pipeline no traslocation .Rdata
# load("../../workfiles/DEG/RNA-seq_DEG_coMMpass_env_noT.Rdata")


# 6.0____________________ PAPER PLOTS ____________________

# Volcano plot of strict DEG 

df <- data.frame(log2_FoldChange=tfit$coefficients[,1],
                 pvalue=tfit$p.value[,1],
                 neg_log10_pvalue= -log10(tfit$p.value[,1]),
                 gene=tfit$genes[,1],
                 strict=dt[,1] %>% as.character())

df$status <- ifelse(df$pvalue<0.05 & df$log2_FoldChange>1, "upreg",
                    ifelse(df$pvalue<0.05 & df$log2_FoldChange < -1, "downreg",NA))

df$status %>% table

df$top20 <- ifelse(df$neg_log10_pvalue> 10, 1,0)
df$top20 %>% table

df$neg_log10_pvalue %>% sort(decreasing = T)

library(ggrepel)
min_log_qval <- df %>% filter(strict!=0) %>% .$neg_log10_pvalue %>% min
df %>% ggplot(aes(x = log2_FoldChange, y = neg_log10_pvalue)) + geom_point(aes(colour=status)) + 
  geom_point(data= df %>% filter(strict!=0), shape=1, colour="black") +
  geom_hline(yintercept = -log10(0.05),linetype=2) +
  geom_hline(yintercept = min_log_qval, linetype=3) +
  geom_vline(xintercept = c(-1,1), linetype=2) +
  ggrepel::geom_text_repel(data = df %>% filter(top20==1), aes(label=gene), max.overlaps = 30)+
  theme_minimal() +
  xlab("Log2 Fold-change")+
  ylab("-log10(p-value)")


dir.create("../../plots/DEG_analysis/", recursive = T, showWarnings = F)
ggsave("../../plots/DEG_analysis/Volcano_plot_strict_DEG.pdf", device = "pdf", width = 8, height = 7)

# no trasnlocations
# ggsave("../../plots/DEG_analysis/Volcano_plot_strict_DEG_noT.pdf", device = "pdf", width = 8, height = 7)


