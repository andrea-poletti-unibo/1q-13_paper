library(tidyverse)

ai_bolo <- data.table::fread("workfiles/AncestralityIndex_table_Bologna_060821.txt")
ai_comm <- data.table::fread("workfiles/AncestralityIndex_table_CoMMpass_0600821.txt")

ai_bolo$alteration == ai_comm$alteration    

# ai_comm$alteration <- recode(ai_comm$alteration, 
#        SeqWGS_WHSC1_CALL ="t_4_14", 
#        SeqWGS_CCND1_CALL="t_11_14",
#        SeqWGS_CCND3_CALL="t_6_14",
#        SeqWGS_MAF_CALL="t_14_16",
#        SeqWGS_MAFB_CALL="t_14_20",)

# ai_comm <- ai_comm %>% filter(!grepl("Seq",ai_comm$alteration))

merge <- left_join(ai_bolo, ai_comm, by = c("alteration"), suffix = c(".Bolo", ".CoMM"),) %>% 
  arrange(desc(Ancestrality_Index.Bolo))


merge$alteration <- merge$alteration %>% str_remove("maj_call_|maj_broad_") %>% str_replace_all("_", " ") %>% str_remove("chr ") %>% str_replace_all("DEL", "Del") %>% str_replace_all("AMP", "Amp")

merge$alteration <- factor(merge$alteration, levels = merge$alteration, ordered = T)

merge$alteration[order(merge$Ancestrality_Index.Bolo)]


gg1 <- merge %>% ggplot(aes(alteration, Ancestrality_Index.Bolo)) +
  geom_point(aes(color="MM-BO"), position = position_nudge(x = -0.15), size=2) +
  geom_point(aes(x=alteration,y=Ancestrality_Index.CoMM, color="CoMMpass"), size=2, position = position_nudge(x = 0.15)) +
  geom_segment(aes(x=alteration, xend=alteration, y=0, yend=Ancestrality_Index.Bolo, color="MM-BO"), position = position_nudge(x = -0.15))+
  geom_segment(aes(x=alteration, xend=alteration, y=0, yend=Ancestrality_Index.CoMM, color="CoMMpass"), position = position_nudge(x = 0.15))+
  ylab("Ancestrality Index") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = c(0.9, 0.8))+
  guides(color=guide_legend(title="Dataset"))

gg1


gg2 <- merge %>% ggplot(aes(Ancestrality_Index.Bolo, Ancestrality_Index.CoMM)) + 
  theme_light() +
  geom_smooth(method = "lm", alpha= 0.5) +
  geom_point(fill="black", colour="black", shape=21, size=3, alpha=0.5) + 
  ggpubr::stat_cor() +
  xlab(label = "Ancestrality Index MM-BO") +
  ylab("Ancestrality Index CoMMpass") +
  ggrepel::geom_text_repel(data = merge %>% filter(Ancestrality_Index.Bolo>2.6 & Ancestrality_Index.CoMM >2.6),
                            aes(label=alteration), force = 1, min.segment.length=0, max.overlaps = 100) +
    theme(legend.position="none", 
          panel.border = element_rect(color = "black",
                                      fill = NA,
                                      size = 1.5))
  
gg2



gg1 + annotation_custom(ggplotGrob(gg2), xmin = 30, xmax = 60, 
                       ymin = 3.1, ymax = 7.5)



ggsave("C:/Users/andre/Alma Mater Studiorum Università di Bologna/PROJECT 1q & 13 - Documenti/Supplementary/Ancestrality_index_defPLOT_231222.pdf", 
       device = "pdf", width = 12, height = 8)
