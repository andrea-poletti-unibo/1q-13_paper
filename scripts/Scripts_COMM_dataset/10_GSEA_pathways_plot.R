
library(tidyverse)

# how to interprert
# https://www.gsea-msigdb.org/gsea/doc/GSEAUserGuideTEXT.htm#_False_Discovery_Rate 


neg <- data.table::fread("input_data/GSEA_results/neg_regulated_pathways.txt")
pos <- data.table::fread("input_data/GSEA_results/pos_regulated_pathways.txt")

neg$class <- "negatively deregulated"
pos$class <- "positively deregulated"


all <- rbind(neg, pos)

RANK <- order(all$NES)

all$Pathway <- all$NAME %>% factor(levels = all$NAME[RANK], ordered = T)

all %>% ggplot(aes(NES, Pathway, colour=class)) + 
  geom_point(aes(size=SIZE, alpha =1-`FDR q-val`)) + 
  facet_wrap(~class, scales = "free_x")+
  scale_alpha_continuous(range = c(0.4,1)) +
  scale_size_continuous(range = c(1,4)) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 6)) 


dir.create("plots/GSEA/", showWarnings = F, recursive = T)

ggsave("plots/GSEA/All_GSEA_pathways.pdf", width = 18, height = 14 )


# selected Pathways

all %>% filter(NAME %>% str_detect("CELL_CYCLE|IMMUNE")) %>% 
  ggplot(aes(NES, Pathway, colour=class)) + 
                 geom_point(aes(size=SIZE, alpha =1-`FDR q-val`)) + 
                 facet_wrap(~class, scales = "free_x") +
  scale_alpha_continuous(range = c(0.3,1)) +
  scale_size_continuous(range = c(2,4)) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 9))

ggsave("plots/GSEA/Selected_GSEA_pathways.pdf", width = 16, height = 8 )
