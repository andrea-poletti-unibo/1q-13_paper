library(tidyverse)
library(data.table)
library(pheatmap)

# load DEG pipeline .Rdata
load("C:/Users/mm_gr/Alma Mater Studiorum Università di Bologna/PROJECT 1q & 13 - Documenti/other/DEG_analysis/7_RNA-seq_DEG_coMMpass_env.Rdata")

# load DEG pipeline no traslocation .Rdata
load("C:/Users/mm_gr/Alma Mater Studiorum Università di Bologna/PROJECT 1q & 13 - Documenti/other/DEG_analysis/7_RNA-seq_DEG_coMMpass_env_noT.Rdata")

# get pts code
vplot<-v
PTS <- colnames(vplot$E)

# import compass data
import <- fread("C:/Users/mm_gr/Alma Mater Studiorum Università di Bologna/PROJECT 1q & 13 - Documenti/complete_database_1q_13_CoMMpass_240121.txt", na.strings=c("","NA"))
import$Study_Visit_iD

# define the MM risk label from commpass broad + focal CNAs
import$AMP_1q_genes_all <- with(import, 
                                ifelse(`AMP_maj_focal_ANP32E`==1 | `AMP_maj_focal_MCL1` ==1 | `AMP_maj_focal_CKS1B`==1, 1,0 ))

import$MMrisk_1q_all <- ifelse( import$`AMP_maj_broad_chr_1q` == 1 | import$AMP_1q_genes_all ==1, 1, 0)
import$MMrisk_13_all <- ifelse( import$`DEL_maj_broad_chr_13q` ==1 | import$`DEL_maj_focal_RB1` ==1, 1, 0)

import$MMrisk_CLASSIFIER_ALL <- 3  - import$MMrisk_1q_all - import$MMrisk_13_all
import$MMrisk_CLASSIFIER_ALL <- dplyr::recode(as.character(import$MMrisk_CLASSIFIER_ALL), "1"="1q&13+","2"="1q/13","3"="1q&13-")

table(import$MMrisk_CLASSIFIER_ALL)

# select only pts with GEP
import_pts <- import %>% filter(Study_Visit_iD %in% PTS)
table(import_pts$MMrisk_CLASSIFIER_ALL)




# define genes

GENES <-c("CCND1","CCND2","CCND3","NSD2","FGFR3", "MAF", "MAFB")
# GENES <-c("CCND1","CCND2")
# GENES <-c("CCND1","CCND2","CCND3")

i <- which(vplot$genes$genes.hgnc_symbol %in% GENES)



# pheatmap(vplot$E[i,], 
#          clustering_method = "ward.D2",
#          scale = "row", 
#          cluster_rows = F,
#          labels_row = vplot$genes$genes.hgnc_symbol[i])


#_____ generate CCND1/2/3 clust label ______ 

GENES <-c("CCND1","CCND2","CCND3")
i <- which(vplot$genes$genes.hgnc_symbol %in% GENES)

mat_c <- vplot$E[i,]
rownames(mat_c) <- vplot$genes$genes.hgnc_symbol[i]

ph <- pheatmap(mat_c, 
               clustering_method = "ward.D2",
               scale = "row", 
               cluster_rows = F)


CLUST <- ph$tree_col %>% cutree(k = 3)
CLUST[CLUST==1] <- "CCND1"
CLUST[CLUST==2] <- "CCND2"
CLUST[CLUST==3] <- "CCND3"

anno <- data.frame(
  cyclin_cluster=CLUST %>% as.factor()
)
rownames(anno) <- PTS

ph <- pheatmap(mat_c, 
               clustering_method = "ward.D2",
               scale = "row", 
               cluster_rows = F,
               annotation_col = anno)




#_______ DEF HEATMAP IN PAPER _______

GENES <-c("CCND1","CCND2","CCND3") # no traslocations

GENES <-c("CCND1","CCND2","CCND3","NSD2","FGFR3", "MAF", "MAFB")

i <- which(vplot$genes$genes.hgnc_symbol %in% GENES)

mat <- vplot$E[i,]
rownames(mat) <- vplot$genes$genes.hgnc_symbol[i]

import_pts$AMP_chr_11 <- ifelse(import_pts$AMP_maj_broad_chr_11p == 1 | import_pts$AMP_all_broad_chr_11q==1, 1,0)

import_pts$AMP_chr_11 %>% table

anno <- data.frame(
  group_1q_13= import_pts$MMrisk_CLASSIFIER_ALL %>% as.factor(),
  cyclin_cluster=CLUST %>% as.factor(),
  t_4_14=import_pts$SeqWGS_WHSC1_CALL %>% as.factor(),
  t_11_14=import_pts$SeqWGS_CCND1_CALL %>% as.factor(),
  t_6_14=import_pts$SeqWGS_CCND3_CALL %>% as.factor(),
  t_14_16=import_pts$SeqWGS_MAF_CALL %>% as.factor(),
  t_14_20=import_pts$SeqWGS_MAFB_CALL %>% as.factor(),
  amp1q=import_pts$MMrisk_1q_all %>% as.factor(),
  del13=import_pts$MMrisk_13_all %>% as.factor(),
  # HD= import_pts$HyperDiploidy %>% as.factor(),
  HD_chr11q= import_pts$AMP_maj_broad_chr_11q %>% as.factor()
)
rownames(anno) <- PTS



anno_colors = list(
  group_1q_13 = c(`1q&13+` = "orangered", `1q&13-` = "turquoise3", `1q/13` = "grey60"),
  cyclin_cluster =  c(`CCND1` = "#66c2a5", `CCND2` = "#fc8d62", `CCND3`="#8da0cb"),
  t_4_14 = c(`1` = "purple", `0` = "grey80"),
  t_11_14= c(`1` = "purple", `0` = "grey80"),
  t_6_14= c(`1` = "purple", `0` = "grey80"),
  t_14_16= c(`1` = "purple", `0` = "grey80"),
  t_14_20= c(`1` = "purple", `0` = "grey80"),
  amp1q=   c(`1` = "green", `0` = "grey80"),
  del13=   c(`1` = "red", `0` = "grey80"),
  # HD= c(`1` = "blue", `0` = "grey80"),
  HD_chr11q= c(`1` = "blue", `0` = "grey80")
)


pheatmap(mat, 
         clustering_method = "ward.D2",
         scale = "row", 
         cluster_rows = F,
         annotation_col = anno,
         annotation_colors = anno_colors, 
         show_colnames = F,
         cellheight = 40,
         width = 16,
         height = 8
)



# IMPORTANT STATS and frequencies in PAPER!
# t involving CCND2
import_pts$t_CCND2 <- ifelse(import_pts$SeqWGS_WHSC1_CALL == 1 | import_pts$SeqWGS_MAF_CALL ==1 | import_pts$SeqWGS_MAFB_CALL ==1 , 1 ,0)

import_pts$t_CCND2 %>% table
table(import_pts$MMrisk_CLASSIFIER_ALL, import_pts$t_CCND2) # percentqage of 1q&13+ that carries a t involving CCND2


# event involving CCND1
import_pts$event_CCND1 <- ifelse(import_pts$SeqWGS_CCND1_CALL == 1 | import_pts$AMP_maj_broad_chr_11q ==1 , 1 ,0)

import_pts$event_CCND1 %>% table
table(import_pts$MMrisk_CLASSIFIER_ALL, import_pts$event_CCND1) # percentqage of 1q&13+ that carries a t involving CCND2





# # SAVE
# pheatmap(mat, 
#          clustering_method = "ward.D2",
#          scale = "row", 
#          cluster_rows = F,
#          annotation_col = anno,
#          annotation_colors = anno_colors, 
#          show_colnames = F,
#          cellheight = 40,
#          filename = "C:/Users/andre/Alma Mater Studiorum Università di Bologna/PROJECT 1q & 13 - Documenti/PAPER_FIGURES_OFFICIAL/DEG_analysis/annotated_genes_heatmap.pdf",
#          width = 16,
#          height = 8
#          )


GENES <-c("CCND1","CCND2","CCND3") # no traslocations

GENES <-c("CCND1","CCND2","CCND3","NSD2","FGFR3", "MAF", "MAFB")

i <- which(vplot$genes$genes.hgnc_symbol %in% GENES)

mat <- vplot$E[i,]
rownames(mat) <- vplot$genes$genes.hgnc_symbol[i]



anno <- data.frame(
  group_1q_13= import_pts$MMrisk_CLASSIFIER_ALL %>% as.factor(),
  cyclin_cluster=CLUST %>% as.factor(),
  t_4_14=import_pts$SeqWGS_WHSC1_CALL %>% as.factor(),
  t_11_14=import_pts$SeqWGS_CCND1_CALL %>% as.factor(),
  t_6_14=import_pts$SeqWGS_CCND3_CALL %>% as.factor(),
  t_14_16=import_pts$SeqWGS_MAF_CALL %>% as.factor(),
  t_14_20=import_pts$SeqWGS_MAFB_CALL %>% as.factor(),
  amp1q=import_pts$MMrisk_1q_all %>% as.factor(),
  del13=import_pts$MMrisk_13_all %>% as.factor(),
  HD= import_pts$HyperDiploidy %>% as.factor()
)
rownames(anno) <- PTS



anno_colors = list(
  group_1q_13 = c(`1q&13+` = "orangered", `1q&13-` = "turquoise3", `1q/13` = "grey60"),
  cyclin_cluster =  c(`CCND1` = "#66c2a5", `CCND2` = "#fc8d62", `CCND3`="#8da0cb"),
  t_4_14 = c(`1` = "purple", `0` = "grey80"),
  t_11_14= c(`1` = "purple", `0` = "grey80"),
  t_6_14= c(`1` = "purple", `0` = "grey80"),
  t_14_16= c(`1` = "purple", `0` = "grey80"),
  t_14_20= c(`1` = "purple", `0` = "grey80"),
  amp1q=   c(`1` = "green", `0` = "grey80"),
  del13=   c(`1` = "red", `0` = "grey80"),
  HD= c(`1` = "blue", `0` = "grey80")
)


pheatmap(mat, 
         clustering_method = "ward.D2",
         scale = "row", 
         cluster_rows = F,
         annotation_col = anno,
         annotation_colors = anno_colors, 
         show_colnames = F,
         cellheight = 40,
         width = 16,
         height = 8
)


# # SAVE
# pheatmap(mat, 
#          clustering_method = "ward.D2",
#          scale = "row", 
#          cluster_rows = F,
#          annotation_col = anno,
#          annotation_colors = anno_colors, 
#          show_colnames = F,
#          cellheight = 40,
#          filename = "C:/Users/andre/Alma Mater Studiorum Università di Bologna/PROJECT 1q & 13 - Documenti/PAPER_FIGURES_OFFICIAL/DEG_analysis/annotated_genes_heatmap.pdf",
#          width = 16,
#          height = 8
#          )



# # SAVE no traslocations
# pheatmap(mat,
#          clustering_method = "ward.D2",
#          scale = "row",
#          cluster_rows = F,
#          annotation_col = anno,
#          annotation_colors = anno_colors,
#          show_colnames = F,
#          cellheight = 40,
#          filename = "C:/Users/andre/Alma Mater Studiorum Università di Bologna/PROJECT 1q & 13 - Documenti/PAPER_FIGURES_OFFICIAL/DEG_analysis/annotated_genes_heatmap_no_trasloc.pdf",
#          width = 16,
#          height = 8
#          )





#########################  other analysis  ################################


#___ table 1____
table <- gmodels::CrossTable(anno$group_1q_13, anno$cyclin_cluster, prop.chisq = F, chisq= T)
table_df <- table$t %>% as.data.frame.matrix()
table_plot <- table$t %>% as.data.frame
table_plot %>% ggplot(aes(y,Freq, fill=x)) + geom_bar(stat = "identity", position = "dodge")



#___ table 2____
anno$traslocation <- ifelse(anno$t_4_14 ==1, "t_4_14",
                            ifelse(anno$t_11_14 ==1, "t_11_14",
                                   ifelse(anno$t_6_14== 1, "t_6_14", 
                                          ifelse(anno$t_14_16==1,"t_14_16",
                                                 ifelse(anno$t_14_20==1,"t_14_20","no_traslocation")))))
anno$type <- ifelse(anno$HD==1, paste0("HD_", anno$traslocation), anno$traslocation)

table <- gmodels::CrossTable(anno$traslocation, anno$cyclin_cluster, prop.chisq = F, chisq= T)
table_df <- table$t %>% as.data.frame.matrix()
table_plot <- table$t %>% as.data.frame
table_plot %>% ggplot(aes(y,Freq, fill=x)) + geom_bar(stat = "identity", position = "dodge")


#___ table 3____
table <- gmodels::CrossTable(anno$traslocation, anno$group_1q_13, prop.chisq = F, chisq= T)
table_df <- table$t %>% as.data.frame.matrix()
table_plot <- table$t %>% as.data.frame
table_plot %>% ggplot(aes(y,Freq, fill=x)) + geom_bar(stat = "identity", position = "dodge")
table_plot %>% filter(x != "no_traslocation") %>%  ggplot(aes(y,Freq, fill=x)) + geom_bar(stat = "identity", position = "dodge")




# _____ stats _______

anno %>% group_by(group_1q_13, cyclin_cluster, t_4_14) %>% summarise(n())

anno %>% filter(cyclin_cluster=="CCND1") %>% .$group_1q_13 %>% table
anno %>% filter(cyclin_cluster=="CCND2") %>% .$group_1q_13 %>% table


anno %>% filter(group_1q_13=="1q&13+") %>% .$cyclin_cluster %>% table
anno %>% filter(group_1q_13=="1q&13-") %>% .$cyclin_cluster %>% table
anno %>% filter(group_1q_13=="1q/13") %>% .$cyclin_cluster %>% table



anno_df <- anno

library(fastDummies)
anno_df <- fastDummies::dummy_columns(anno_df)

library(nnet)

m <- glm(cyclin_cluster_CCND2~`group_1q_13_1q&13+`, family="binomial", data=anno_df)
summary(m)

m <- glm(cyclin_cluster_CCND1~`group_1q_13_1q&13-`, family="binomial", data=anno_df)
summary(m)

m <- glm(cyclin_cluster_CCND3~`group_1q_13_1q/13`, family="binomial", data=anno_df)
summary(m)


library(kableExtra)

table_df %>% kbl(caption = "1q&13 groups and Cyclin D clusters") %>%
  kable_classic(full_width = F, html_font = "Cambria") %>% 
  footnote( general_title = "Pearson χ² test p < 0.0001", general = "", footnote_as_chunk = T) %>% 
  save_kable( "C:/Users/mm_gr/Alma Mater Studiorum Università di Bologna/PROJECT 1q & 13 - Documenti/PAPER_FIGURES_OFFICIAL/DEG_analysis/table_MMrisk_CyclinClusters.pdf",vwidth=500, vheight=300)




