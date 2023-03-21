
library(tidyverse)

# tab_down <- readxl::read_xlsx("C:/Users/andre/Alma Mater Studiorum Università di Bologna/PROJECT 1q & 13 - Documenti/SUBMISSION_Cancer_Discovery/SUPPLEMENTARY_ITEMS/Sup_table_S4.xlsx", sheet = 1)
# 
# tab_up <- readxl::read_xlsx("C:/Users/andre/Alma Mater Studiorum Università di Bologna/PROJECT 1q & 13 - Documenti/SUBMISSION_Cancer_Discovery/SUPPLEMENTARY_ITEMS/Sup_table_S4.xlsx", sheet = 2)
# 
# dir.create("data/GSEA_official_res")
# 
# write_tsv(tab_down,"data/GSEA_official_res/neg_regulated_pathways.txt")
# write_tsv(tab_up,"data/GSEA_official_res/pos_regulated_pathways.txt")


# how to interprert
# https://www.gsea-msigdb.org/gsea/doc/GSEAUserGuideTEXT.htm#_False_Discovery_Rate 


neg <- data.table::fread("data/GSEA_official_res/neg_regulated_pathways.txt")
pos <- data.table::fread("data/GSEA_official_res/pos_regulated_pathways.txt")

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


# dir.create("C:/Users/andre/Alma Mater Studiorum Università di Bologna/PROJECT 1q & 13 - Documenti/PAPER_FIGURES_OFFICIAL/GSEA_pathways")

ggsave("C:/Users/andre/Alma Mater Studiorum Università di Bologna/PROJECT 1q & 13 - Documenti/PAPER_FIGURES_OFFICIAL/GSEA_pathways/All_GSEA_pathways.pdf", width = 18, height = 14 )


# selected Pathways

all %>% filter(NAME %>% str_detect("CELL_CYCLE|IMMUNE")) %>% 
  ggplot(aes(NES, Pathway, colour=class)) + 
                 geom_point(aes(size=SIZE, alpha =1-`FDR q-val`)) + 
                 facet_wrap(~class, scales = "free_x") +
  scale_alpha_continuous(range = c(0.3,1)) +
  scale_size_continuous(range = c(2,4)) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 9))

ggsave("C:/Users/andre/Alma Mater Studiorum Università di Bologna/PROJECT 1q & 13 - Documenti/PAPER_FIGURES_OFFICIAL/GSEA_pathways/Selected_GSEA_pathways.pdf", width = 16, height = 8 )
