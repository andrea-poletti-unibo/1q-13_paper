library(tidyverse)
library(data.table)

import <- fread("D:/analisi_in_corso/risk_1q_13/CoMMpass/complete_database_1q_13_CoMMpass.txt")

mutations <- fread("D:/CoMMpass_data/IA13a/Somatic_Mutations-SNV_and_INDEL/MMRF_CoMMpass_IA13a_All_Canonical_NS_Variants.txt")

mutations_counts <- fread("D:/CoMMpass_data/IA13a/Somatic_Mutations-SNV_and_INDEL/MMRF_CoMMpass_IA13a_All_Canonical_NS_Variants_ENSG_Mutation_Counts.txt")

# define the MM risk label from commpass broad + focal CNAs
import$AMP_1q_genes_all <- with(import, 
                                ifelse(`AMP_maj_focal_ANP32E`==1 | `AMP_maj_focal_MCL1` ==1 | `AMP_maj_focal_CKS1B`==1, 1,0 ))

import$MMrisk_1q_all <- ifelse( import$`AMP_maj_broad_chr_1q` == 1 | import$AMP_1q_genes_all ==1, 1, 0)
import$MMrisk_13_all <- ifelse( import$`DEL_maj_broad_chr_13q` ==1 | import$`DEL_maj_focal_RB1` ==1, 1, 0)

import$MMrisk_CLASSIFIER_ALL <- 3  - import$MMrisk_1q_all - import$MMrisk_13_all
import$MMrisk_CLASSIFIER_ALL <- dplyr::recode(as.character(import$MMrisk_CLASSIFIER_ALL), "1"="risk1","2"="risk2","3"="risk3")

table(import$MMrisk_CLASSIFIER_ALL)


mut_count_mat <- mutations_counts[,-1] %>% as.matrix()

mut_burden<- apply(mut_count_mat, 2 , sum)

mut_burden

mut_burden_per_pz <- data.frame(Study_Visit_iD= colnames(mut_count_mat),
                                burden = mut_burden)

mut_burden_per_pz$burden %>% sort %>% plot
mut_burden_per_pz$burden %>% summary


merge <- left_join(import, mut_burden_per_pz, by="Study_Visit_iD")


res <- merge %>% group_by(MMrisk_CLASSIFIER_ALL) %>% summarise(med_burden= median(burden, na.rm = T),
                                                        mean_burden= mean(burden, na.rm = T),
                                                        sd_burden=sd(burden, na.rm = T),
                                                        n_NAs= is.na(burden) %>% sum)

res[,1:3] %>% View

burdens_risk1 <- merge$burden[merge$MMrisk_CLASSIFIER_ALL=="risk1"]
burdens_risk3 <- merge$burden[merge$MMrisk_CLASSIFIER_ALL=="risk3"]
burdens_risk2and3 <- merge$burden[merge$MMrisk_CLASSIFIER_ALL !="risk1"]
burdens_risk1and2 <- merge$burden[merge$MMrisk_CLASSIFIER_ALL !="risk3"]





t.test(burdens_risk1, burdens_risk2and3) # p-value = 0.01785

t.test(burdens_risk3, burdens_risk1) # p-value = 0.03101

t.test(burdens_risk3, burdens_risk1and2) # p-value = 0.1786






merge %>% ggplot(aes(MMrisk_CLASSIFIER_ALL, burden)) + 
  geom_violin() + 
  geom_boxplot() + 
  ylim(0,250) + 
  geom_jitter(width = 0.1)

