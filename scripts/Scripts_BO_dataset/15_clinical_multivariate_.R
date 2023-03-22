library(data.table)
library(tidyverse)
library(survival)
library(survminer)
library(RODBC)
library(broom)
options(scipen=999)


# ======== import clinical data already created by 230321 ======== 
import <- data.table::fread("results/complete_database_BO_1q_13.txt") %>% as.data.frame()
names(import)

outpath <- "results/Clinical_anlysis/BO_dataset/Multivariate/"
dir.create(outpath)

#========== MM RISK CREATION ============

import$AMP_1q_genes_all <- with(import, 
                                ifelse(`AMP_maj-focal_ANP32E`==1 | `AMP_maj-focal_MCL1` ==1 | `AMP_maj-focal_CKS1B`==1, 1,0 ))

import$MMrisk_1q_all<- ifelse( import$`AMP_maj-broad_chr_1q` == 1 | import$AMP_1q_genes_all ==1, 1, 0)
import$MMrisk_13_all <- ifelse( import$`DEL_maj-broad_chr_13q` ==1 | import$`DEL_maj-focal_RB1` ==1, 1, 0)

import$MMrisk_CLASS <- 3  - import$MMrisk_1q_all - import$MMrisk_13_all

import$MMrisk_AB_ALL <- ifelse(import$MMrisk_CLASS == 2 & import$MMrisk_1q_all ==1,
                               "2_1q", 
                               ifelse(import$MMrisk_CLASS == 2 & import$MMrisk_13_all ==1,
                                      "2_13",
                                      import$MMrisk_CLASS))

import$MMrisk_class <- paste0("risk_",import$MMrisk_CLASS) %>% as.factor() %>% relevel("risk_3")
import$MMrisk_class <- import$MMrisk_class %>% recode_factor(risk_1 ="1q&13+", risk_2 ="1q/13", risk_3 ="1q&13-") %>% relevel("1q&13-")

import$MMrisk_allclass <- ifelse(import$MMrisk_class=="1q/13" & import$MMrisk_1q_all==1, "gain 1q only ", 
                                 ifelse(import$MMrisk_class=="1q/13" & import$MMrisk_13_all==1, "del 13q only",import$MMrisk_class %>% as.character())) %>% as.factor()

import$MMrisk_allclass


import$MMrisk1 <- ifelse(import$MMrisk_CLASS==1, 1, 0)
import$MMrisk2 <- ifelse(import$MMrisk_CLASS==2, 1, 0)
import$MMrisk3 <- ifelse(import$MMrisk_CLASS==3, 1, 0)

import$risk1_VS_1q_VS_other <- recode(import$MMrisk_AB_ALL, "1"="1q&13", "2_1q"="1q_only", "2_13"="other", "3"="other")

#___ ultra MMrisk ____
import$ULTRA_MMrisk <- ifelse(import$ISS == 3 & import$MMrisk_CLASS == 1,
                              "Ultra_High",
                              ifelse(import$ISS==1 & import$MMrisk_CLASS == 3,
                                     "Ultra_Low",
                                     "Other")) %>% as.factor()
levels(import$ULTRA_MMrisk)
import$ULTRA_MMrisk %>% table

import$ULTRA_High <- ifelse(import$ULTRA_MMrisk=="Ultra_High",1,0)
import$ULTRA_Low  <- ifelse(import$ULTRA_MMrisk=="Ultra_Low",1,0)

#_____ create pure MMrisk1 and pure CCND2 vars _____

import$CCND2_traslocation <- 0
import$CCND2_traslocation[import$FISH_T_4_14==1] <- 1
import$CCND2_traslocation[import$FISH_T_14_16==1] <- 1
import$CCND2_traslocation[import$FISH_T_14_20==1] <- 1
import$CCND2_traslocation[is.na(import$FISH_T_4_14) & is.na(import$FISH_T_14_16) & is.na(import$FISH_T_14_20)] <- NA


import$MMrisk_class_t_CCND2 <- ifelse(import$MMrisk_CLASS==1 & import$CCND2_traslocation ==1 , "t&1q&13+",
                                      ifelse(import$MMrisk_CLASS==1, "1q&13+_pure", import$MMrisk_allclass %>% as.character))


import$FISH_T_4_14_pure <- ifelse(import$MMrisk_CLASS != 1 & import$FISH_T_4_14 ==1, 1, 0)
import$FISH_T_14_16_pure <- ifelse(import$MMrisk_CLASS != 1 & import$FISH_T_14_16 ==1, 1, 0)
import$FISH_T_14_20_pure <- ifelse(import$MMrisk_CLASS != 1 & import$FISH_T_14_20 ==1, 1, 0)




#___ ultra MMrisk CCND2 Traslocations ____
import$CCND2_ULTRA <- ifelse(import$R_ISS == 3 & import$MMrisk_class_t_CCND2 == "t&1q&13+",
                              "Ultra_t_High",
                              ifelse(import$R_ISS==1 & import$MMrisk_class_t_CCND2 == "1q&13-",
                                     "Ultra_Low",
                                     "Other")) %>% as.factor()
levels(import$CCND2_ULTRA)
import$CCND2_ULTRA %>% table


#========== create additional variables ========== 

import$ALBUMIN_m_3.5 <- ifelse(import$ALBUMIN < 3.5, 1,0)

import$B2M_level <- ifelse(import$B2M < 3.5, "<3.5", ifelse(import$B2M > 5.5, ">5.5", "normal")) %>% as.factor()
import$B2M_level <- import$B2M_level %>% relevel(ref = 3)
levels(import$B2M_level)


import$Induction_therapy <- ifelse(!import$INDUCTION_THERAPY_FRONT_LINE %in% c("VCD", "VTD", "TD" ), "other", import$INDUCTION_THERAPY_FRONT_LINE %>% as.character() )
import$Induction_therapy %>% table

import$Induction_response_M_VGPR <- ifelse(import$INDUCTION_RESPONSE %in% c("VGPR", "nCR", "CR", "sCR"), 1,0)
import$Induction_response_M_CR <- ifelse(import$INDUCTION_RESPONSE %in% c("nCR", "CR", "sCR"), 1,0)

import$ASCT<-import$TX_I_LINE_1_0
import$MAINTENANCE<-import$MAINTENANCE_YES_NO
import$CONSOLIDATION<-import$CONSOLIDATION_YES_NO

table(import$NOTE_TRANSLOCATION)
import$IgH_translocation_type<-ifelse(import$NOTE_TRANSLOCATION=="non traslocato" ,"no translocation",
                                        ifelse(import$NOTE_TRANSLOCATION=="pannello incompleto", NA, as.character(import$NOTE_TRANSLOCATION))) %>% as.factor()
levels(import$IgH_translocation_type)

#========= additional data management =========

class(import$PROTOCOL_REV)
import$PROTOCOL_REV<-as.factor(import$PROTOCOL_REV)
levels(import$PROTOCOL_REV)

import$OLD <-import$OLD_0_1

class(import$IgH_translocation_type)
levels(import$IgH_translocation_type)

import$MMrisk_AB_ALL<-as.factor(import$MMrisk_AB_ALL)
levels(import$MMrisk_AB_ALL)

class(import$ISS)
import$ISS<-as.factor(import$ISS)
levels(import$ISS)

import$ISS<-relevel(import$ISS, ref=2)
levels(import$ISS)

import$R_ISS <- as.factor(import$R_ISS) %>% relevel(ref = 2)
levels(import$R_ISS)

# ========= create surv =========

OS <- Surv(import$OS_MONTHS, import$OS_EVENT)
PFS <- Surv(import$PFS_I_MONTHS, import$PFS_I_EVENT)



#.................# ggforest2 function def #......................

ggforest2 <- function (model, data = NULL, main = "Hazard ratio", cpositions = c(0.02, 
                                                                                 0.22, 0.4), fontsize = 0.7, refLabel = "reference", noDigits = 2) 
{
  conf.high <- conf.low <- estimate <- NULL
  stopifnot(inherits(model, "coxph"))
  data <- survminer:::.get_data(model, data = data)
  terms <- attr(model$terms, "dataClasses")[-1]
  coef <- as.data.frame(tidy(model, conf.int = TRUE))
  gmodel <- glance(model)
  allTerms <- lapply(seq_along(terms), function(i) {
    var <- names(terms)[i]
    if (var %in% colnames(data)) {
      if (terms[i] %in% c("factor", "character")) {
        adf <- as.data.frame(table(data[, var]))
        cbind(var = var, adf, pos = 1:nrow(adf))
      }
      else if (terms[i] == "numeric") {
        data.frame(var = var, Var1 = "", Freq = nrow(data), 
                   pos = 1)
      }
      else {
        vars = grep(paste0("^", var, "*."), coef$term, 
                    value = TRUE)
        data.frame(var = vars, Var1 = "", Freq = nrow(data), 
                   pos = seq_along(vars))
      }
    }
    else {
      message(var, "is not found in data columns, and will be skipped.")
    }
  })
  allTermsDF <- do.call(rbind, allTerms)
  colnames(allTermsDF) <- c("var", "level", "N", "pos")
  inds <- apply(allTermsDF[, 1:2], 1, paste0, collapse = "")
  rownames(coef) <- gsub(coef$term, pattern = "`", replacement = "")
  toShow <- cbind(allTermsDF, coef[inds, ])[, c("var", "level", 
                                                "N", "p.value", "estimate", "conf.low", "conf.high", 
                                                "pos")]
  toShowExp <- toShow[, 5:7]
  toShowExp[is.na(toShowExp)] <- 0
  toShowExp <- format(exp(toShowExp), digits = noDigits)
  toShowExpClean <- data.frame(toShow, pvalue = signif(toShow[, 
                                                              4], noDigits + 1), toShowExp)
  toShowExpClean$stars <- paste0(round(toShowExpClean$p.value, 
                                       noDigits + 1), " ", ifelse(toShowExpClean$p.value < 
                                                                    0.05, "*", ""), ifelse(toShowExpClean$p.value < 0.01, 
                                                                                           "*", ""), ifelse(toShowExpClean$p.value < 0.001, "*", 
                                                                                                            ""))
  toShowExpClean$ci <- paste0("(", toShowExpClean[, "conf.low.1"], 
                              " - ", toShowExpClean[, "conf.high.1"], ")")
  toShowExpClean$estimate.1[is.na(toShowExpClean$estimate)] = refLabel
  toShowExpClean$stars[which(toShowExpClean$p.value < 0.001)] = "<0.001 ***"
  toShowExpClean$stars[is.na(toShowExpClean$estimate)] = ""
  toShowExpClean$ci[is.na(toShowExpClean$estimate)] = ""
  toShowExpClean$estimate[is.na(toShowExpClean$estimate)] = 0
  toShowExpClean$var = as.character(toShowExpClean$var)
  toShowExpClean$var[duplicated(toShowExpClean$var)] = ""
  toShowExpClean$N <- paste0("(N=", toShowExpClean$N, ")")
  toShowExpClean <- toShowExpClean[nrow(toShowExpClean):1, 
  ]
  rangeb <- range(toShowExpClean$conf.low, toShowExpClean$conf.high, 
                  na.rm = TRUE)
  breaks <- axisTicks(rangeb/2, log = TRUE, nint = 7)
  rangeplot <- rangeb
  rangeplot[1] <- rangeplot[1] - diff(rangeb)
  rangeplot[2] <- rangeplot[2] + 0.15 * diff(rangeb)
  width <- diff(rangeplot)
  y_variable <- rangeplot[1] + cpositions[1] * width
  y_nlevel <- rangeplot[1] + cpositions[2] * width
  y_cistring <- rangeplot[1] + cpositions[3] * width
  y_stars <- rangeb[2]
  x_annotate <- seq_len(nrow(toShowExpClean))
  annot_size_mm <- fontsize * as.numeric(grid::convertX(unit(theme_get()$text$size, 
                                                             "pt"), "mm"))
  p <- ggplot(toShowExpClean, aes(seq_along(var), exp(estimate))) + 
    geom_rect(aes(xmin = seq_along(var) - 0.5, xmax = seq_along(var) + 
                    0.5, ymin = exp(rangeplot[1]), ymax = exp(rangeplot[2]), 
                  fill = ordered(seq_along(var)%%2 + 1))) + scale_fill_manual(values = c("#FFFFFF33", 
                                                                                         "#00000033"), guide = "none") + geom_point(pch = 15, 
                                                                                                                                    size = 4) + geom_errorbar(aes(ymin = exp(conf.low), 
                                                                                                                                                                  ymax = exp(conf.high)), width = 0.15) + geom_hline(yintercept = 1, 
                                                                                                                                                                                                                     linetype = 3) + coord_flip(ylim = exp(rangeplot)) + 
    ggtitle(main) + scale_y_log10(name = "", labels = sprintf("%g", 
                                                              breaks), expand = c(0.02, 0.02), breaks = breaks) + 
    theme_light() + theme(panel.grid.minor.y = element_blank(), 
                          panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), 
                          legend.position = "none", panel.border = element_blank(), 
                          axis.title.y = element_blank(), axis.text.y = element_blank(), 
                          axis.ticks.y = element_blank(), plot.title = element_text(hjust = 0.5)) + 
    xlab("") + annotate(geom = "text", x = x_annotate, y = exp(y_variable), 
                        label = toShowExpClean$var, fontface = "bold", hjust = 0, 
                        size = annot_size_mm) + annotate(geom = "text", x = x_annotate, 
                                                         y = exp(y_nlevel), hjust = 0, label = toShowExpClean$level, 
                                                         vjust = -0.1, size = annot_size_mm) + annotate(geom = "text", 
                                                                                                        x = x_annotate, y = exp(y_nlevel), label = toShowExpClean$N, 
                                                                                                        fontface = "italic", hjust = 0, vjust = ifelse(toShowExpClean$level == 
                                                                                                                                                         "", 0.5, 1.1), size = annot_size_mm) + annotate(geom = "text", 
                                                                                                                                                                                                         x = x_annotate, y = exp(y_cistring), label = toShowExpClean$estimate.1, 
                                                                                                                                                                                                         size = annot_size_mm, vjust = ifelse(toShowExpClean$estimate.1 == 
                                                                                                                                                                                                                                                "reference", 0.5, -0.1)) + annotate(geom = "text", 
                                                                                                                                                                                                                                                                                    x = x_annotate, y = exp(y_cistring), label = toShowExpClean$ci, 
                                                                                                                                                                                                                                                                                    size = annot_size_mm, vjust = 1.1, fontface = "italic") + 
    annotate(geom = "text", x = x_annotate, y = exp(y_stars), 
             label = toShowExpClean$stars, size = annot_size_mm, 
             hjust = -0.2, fontface = "italic") + annotate(geom = "text", 
                                                           x = 0.5, y = exp(y_variable), label = paste0("# Events: ", 
                                                                                                        gmodel$nevent, "; Global p-value (Log-Rank): ", 
                                                                                                        format.pval(gmodel$p.value.log, eps = ".001"), " \nAIC: ", 
                                                                                                        round(gmodel$AIC, 2), "; Concordance Index: ", round(gmodel$concordance, 
                                                                                                                                                             2)), size = annot_size_mm, hjust = 0, vjust = 1.2, 
                                                           fontface = "italic")
  gt <- ggplot_gtable(ggplot_build(p))
  gt$layout$clip[gt$layout$name == "panel"] <- "off"
  ggpubr::as_ggplot(gt)
}

#..........................................................


######################### MULTIVARIATE ##############################

#===========  MODEL 1 - multivariate MMrisk =============== 

# OS model 1
mv_R <- coxph(OS ~ OLD + SEX + 
                ISS + 
                HB_m_105 + PLT_m_150 +  
                MMrisk_allclass +
                FISH_Del_17p + FISH_Del_1p36 + HyperDiploidy +
                IgH_translocation_type +
                strata(PROTOCOL_REV), data = import )
mv_R %>% summary()

ggforest2(mv_R, main = "OS MM-BO (model 1)") 

ggsave(paste0(outpath,"FOREST_OS_Bolo_model1.pdf"), width = 10, height = 7)


# PFS model 1
mv_R <- coxph(PFS ~ OLD + SEX + 
                ISS +
                HB_m_105 + PLT_m_150 +  
                MMrisk_allclass +
                FISH_Del_17p + FISH_Del_1p36 + HyperDiploidy +
                IgH_translocation_type +
                strata(PROTOCOL_REV), data = import )
mv_R %>% summary()

ggforest2(mv_R, main = "PFS MM-BO (model 1)") 

ggsave(paste0(outpath,"FOREST_PFS_Bolo_model1.pdf"), width = 10, height = 7)


#=========== MODEL 2 - multivariate MMrisk and CCND2 t =============== 

# OS model 2
mv_R <- coxph(OS ~ OLD + SEX + 
                ISS +
                HB_m_105 + PLT_m_150 +  
                MMrisk_class_t_CCND2 +
                FISH_Del_17p + FISH_Del_1p36 + HyperDiploidy +
                FISH_T_4_14_pure + FISH_T_14_16_pure + FISH_T_14_20_pure + FISH_T_11_14 + FISH_T_6_14 +
                strata(PROTOCOL_REV), data = import )
mv_R %>% summary()

ggforest2(mv_R, main = "OS MM-BO (model 2)") 

ggsave(paste0(outpath,"FOREST_OS_Bolo_model2.pdf"), width = 10, height = 7)


# PFS model 2
mv_R <- coxph(PFS ~ OLD + SEX + 
                ISS +
                HB_m_105 + PLT_m_150 +  
                MMrisk_class_t_CCND2 +
                FISH_Del_17p + FISH_Del_1p36 + HyperDiploidy +
                FISH_T_4_14_pure + FISH_T_14_16_pure + FISH_T_14_20_pure + FISH_T_11_14 + FISH_T_6_14 +
                strata(PROTOCOL_REV), data = import )
mv_R %>% summary()

ggforest2(mv_R, main = "PFS MM-BO (model 2)") 

ggsave(paste0(outpath,"FOREST_PFS_Bolo_model2.pdf"), width = 10, height = 7)



