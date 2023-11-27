library(tidyverse)
library(data.table)
library(ggpubr)


import<-fread("results/complete_database_BO_1q_13.txt")

rownames(import) <- import$SNP_FILE_NAME %>% str_replace(" ", "_")
import <- rownames_to_column(import, "sample")

# select only relevant variables to associate
names(import)
data<-import %>% select(ID, contains("maj_broad"), contains("maj_focal"), contains("maj_call"), HyperDiploidy, matches("^t_"))

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

# exclude calls from HyperDiploid chromosomes and single translocations for MDS dataset creation 
data <- data %>% select( -c(AMP_maj_broad_chr_3p, 
                            AMP_maj_broad_chr_3q, 
                            AMP_maj_broad_chr_5p, 
                            AMP_maj_broad_chr_5q, 
                            AMP_maj_broad_chr_7p, 
                            AMP_maj_broad_chr_7q,
                            AMP_maj_broad_chr_9p, 
                            AMP_maj_broad_chr_9q,
                            AMP_maj_broad_chr_11p, 
                            AMP_maj_broad_chr_11q,
                            AMP_maj_broad_chr_15q,
                            AMP_maj_broad_chr_19p,
                            AMP_maj_broad_chr_19q,
                            AMP_maj_broad_chr_21q,
                            t_4_14,
                            t_11_14,
                            t_others
))

# renaming calls column
names(data) <- names(data) %>% str_replace("maj_broad_chr_|maj_call_chr_","")
names(data)




# filter data for missing data (missing traslocation)


missing_T <- is.na(import$t_11_14) + is.na(import$t_4_14) + is.na(import$t_6_14) + is.na(import$t_14_16) + is.na(import$t_14_20)
table(missing_T)
import$missing_T <- missing_T

df.2<- filter(data, !is.na(t_IgH))

df.2 <- mutate(df.2, SUM=rowSums(df.2[,-1]))


df.2$group <- ifelse( df.2$AMP_1q == 1 & df.2$DEL_13q ==1, "1q&13+",
                      ifelse(df.2$AMP_1q == 1 | df.2$DEL_13q ==1, "1q/13", "1q&13-") ) 

df.2$group %>% table
df.2$SUM

df.2$SUM %>% summary

  
  
KS <- kruskal.test(df.2$group, df.2$SUM)

KS$p.value

GC_plus <- df.2 %>% filter(group=="1q&13+") %>% .$SUM
GC_plus_other <- df.2 %>% filter(group != "1q&13+") %>% .$SUM

df.2$group1 <- ifelse(df.2$group == "1q&13+", "1q&13+","other")
kruskal.test(df.2$group1, df.2$SUM)

df.2$group2 <- ifelse(df.2$group == "1q/13", "1q/13","other")
kruskal.test(df.2$group2, df.2$SUM)

df.2$group3 <- ifelse(df.2$group == "1q&13-", "1q&13-","other")
kruskal.test(df.2$group3, df.2$SUM)



df.2$group <- df.2$group %>% factor(levels = c("1q&13-","1q/13","1q&13+"))

df.2 %>% ggplot(aes(group, SUM, colour=group)) +
  geom_boxplot() +
  geom_jitter(width = 0.2) +
  ggpubr::stat_compare_means()


df.2 %>% ggplot(aes(group2, SUM, colour=group2)) +
  geom_boxplot() +
  geom_jitter(width = 0.2) +
  ggpubr::stat_compare_means()

wilcox.test()
df.2 %>% ggplot(aes(group1, SUM, colour=group1)) +
  geom_boxplot() +
  geom_jitter(width = 0.2) +
  ggpubr::stat_compare_means()

df.2 %>% ggplot(aes(group3, SUM, colour=group3)) +
  geom_boxplot() +
  geom_jitter(width = 0.2) +
  ggpubr::stat_compare_means()


#________ def plots _____ paper


my_comparisons <- list( c("1q&13-", "1q/13"), 
                        c("1q/13", "1q&13+"), 
                        c("1q&13-", "1q&13+") )


gg <- df.2 %>% ggplot(aes(group, SUM)) +
  geom_violin(aes(fill=group, colour=group), alpha=0.5) +   
  geom_boxplot(aes( colour=group), alpha=0.5, width=0.5, outlier.shape = NA) +
  geom_jitter(width = 0.12, alpha=0.3) + 
  stat_summary(geom="text", fun =quantile,
               aes(label=sprintf("%1.2f", after_stat(y)), color=group),
               position=position_nudge(x= 0.5), size=3.5) +
  stat_compare_means(comparisons = my_comparisons)+ 
  stat_compare_means(label.y.npc =0.75, label.x.npc = 0) + 
  xlab("") + ylab("Genomic Complexity") +
  # custom colors
  scale_fill_manual(values=c("turquoise3", "grey60", "orangered")) +
  scale_color_manual(values=c("turquoise3", "grey60", "orangered")) +
  # legend bottom
  theme(legend.position = "bottom") 

gg


dir.create("plots/Genomic_complexity")

write_tsv(df.2, "plots/Genomic_complexity/genomic_complexity_data.txt")

ggsave(filename = "plots/Genomic_complexity/genomic_complexity_in_MMriskGroups.pdf", 
       units = "in", dpi = 300, width = 6, height = 6)

ggsave(filename = "plots/Genomic_complexity/genomic_complexity_in_MMriskGroups.svg", 
       units = "in", dpi = 300, width = 6, height = 6)


df.2 %>% group_by(group) %>% summarize(mean_complexity=mean(SUM), 
                                       CI_lower=SUM %>% gmodels::ci() %>% .[2],
                                       CI_upper=SUM %>% gmodels::ci() %>% .[3])
