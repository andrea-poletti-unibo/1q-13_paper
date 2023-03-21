options(scipen=999)
library(tidyverse)
library(survival)
library(survminer)
library(RODBC)
library(broom)

#______import clinical data already created
import <- data.table::fread("clinical_data_COMMPASS-IA13_210122.txt") %>% as.data.frame()

outpath <- "C:/Users/andre/Alma Mater Studiorum Università di Bologna/PROJECT 1q & 13 - Documenti/PAPER_FIGURES_OFFICIAL/forest_plots/"

#_____ADDITIONAL DATA MANAGEMENT____
import$LDH_level<-ifelse(import$LDH_level=="",NA, import$LDH_level)
import$LDH_level<-as.factor(import$LDH_level)
levels(import$LDH_level)

import$R_ISS<-as.factor(import$R_ISS)
levels(import$R_ISS)

import$SEX <- ifelse(import$D_PT_gender==1, "M", "F") %>% as.factor()

import$SEX<-as.factor(import$SEX)
import$SEX<-relevel(import$SEX, ref = 2)

import$LDH <- import$D_LAB_chem_ldh

import$ASCT <- import$sctflag

import$FIRST_THERAPHY

import$ALBUMIN <- import$D_LAB_chem_albumin
import$ALBUMIN_m_3.5 <- ifelse(import$ALBUMIN < 35, 1,0)

import$B2M <- import$D_LAB_serum_beta2_microglobulin
import$B2M_level <- ifelse(import$B2M < 3.5, "<3.5", ifelse(import$B2M > 5.5, ">5.5", "normal")) %>% as.factor()
import$B2M_level <- import$B2M_level %>% relevel(ref = 3)
levels(import$B2M_level)

import$SeqFISH_Del_TP53 <- import$DEL_maj_focal_TP53

import$SeqFISH_Del_1p_CDKN2C <- ifelse(import$DEL_maj_focal_CDKN2C ==1, 1,0)


#============ create trasnlocations vars ============
import$IgH_translocation_type <-ifelse(import$SeqWGS_WHSC1_CALL==1, "t(4;14)",
                                       ifelse(import$SeqWGS_CCND1_CALL==1, "t(11;14)",
                                              ifelse(import$SeqWGS_CCND3_CALL==1,"t(6;14)",
                                                     ifelse(import$SeqWGS_MAF_CALL==1,"t(14;16)",
                                                            ifelse(import$SeqWGS_MAFB_CALL==1,"t(14;20)", 
                                                                   "no translocation")))))

table(import$IgH_translocation_type)
table(import$IgH_translocation_type)/840

import$IgH_translocation_type<-as.factor(import$IgH_translocation_type)
levels(import$IgH_translocation_type)

import$SeqFISH_T_4_14 <- ifelse(import$SeqWGS_WHSC1_CALL==1, 1, 0)
import$SeqFISH_T_11_14 <- ifelse(import$SeqWGS_CCND1_CALL==1, 1, 0)
import$SeqFISH_T_14_16 <- ifelse(import$SeqWGS_MAF_CALL==1, 1, 0)
import$SeqFISH_T_14_20 <- ifelse(import$SeqWGS_MAFB_CALL==1, 1, 0)
import$SeqFISH_T_6_14 <- ifelse(import$SeqWGS_CCND3_CALL==1, 1, 0)


#========== MM RISK CREATION ============
import$AMP_1q_genes_all <- with(import, 
                                ifelse(`AMP_maj_focal_ANP32E`==1 | `AMP_maj_focal_MCL1` ==1 | `AMP_maj_focal_CKS1B`==1, 1,0 ))

import$MMrisk_1q_all<- ifelse( import$`AMP_maj_broad_chr_1q` == 1 | import$AMP_1q_genes_all ==1, 1, 0)
import$MMrisk_13_all <- ifelse( import$`DEL_maj_broad_chr_13q` ==1 | import$`DEL_maj_focal_RB1` ==1, 1, 0)

import$MMrisk_CLASS <- 3  - import$MMrisk_1q_all - import$MMrisk_13_all
import$MMrisk_class <- paste0("risk_",import$MMrisk_CLASS) %>% as.factor() %>% relevel("risk_3")

import$MMrisk_class <- import$MMrisk_class %>% recode_factor(risk_1 ="1q&13+", risk_2 ="1q/13", risk_3 ="1q&13-") %>% relevel("1q&13-")
import$MMrisk_class %>% table()

import$MMrisk_allclass <- ifelse(import$MMrisk_class=="1q/13" & import$MMrisk_1q_all==1, "gain 1q only ", 
                                 ifelse(import$MMrisk_class=="1q/13" & import$MMrisk_13_all==1, "del 13q only",import$MMrisk_class %>% as.character())) %>% as.factor()

import$MMrisk_allclass

import$MMrisk1 <- ifelse(import$MMrisk_CLASS==1, 1, 0)
import$MMrisk2 <- ifelse(import$MMrisk_CLASS==2, 1, 0)
import$MMrisk3 <- ifelse(import$MMrisk_CLASS==3, 1, 0)

import$MMrisk_AB_ALL <- ifelse(import$MMrisk_CLASS == 2 & import$MMrisk_1q_all ==1,
                               "2_1q", 
                               ifelse(import$MMrisk_CLASS == 2 & import$MMrisk_13_all ==1,
                                      "2_13",
                                      import$MMrisk_CLASS))

import$MMrisk_AB_ALL<-as.factor(import$MMrisk_AB_ALL)
import$MMrisk_AB_ALL<-relevel(import$MMrisk_AB_ALL, ref=4)
levels(import$MMrisk_AB_ALL)

import$risk1_VS_1q_VS_other <- recode(import$MMrisk_AB_ALL, "1"="1q&13", "2_1q"="1q_only", "2_13"="other", "3"="other")

#_____ COMPUTE ISS ______

import$ISS <- ifelse(import$D_LAB_serum_beta2_microglobulin>= 5.5, 3,
                     ifelse(import$D_LAB_serum_beta2_microglobulin< 3.5 & import$D_LAB_chem_albumin >= 35, 1,2)) %>% as.factor() %>% relevel(ref = 2)
import$ISS
table(import$ISS)


import$R_ISS <- import$R_ISS %>% relevel(ref = 2)


#===== create pure vars MMrisk and pure CCND2 =======

import$CCND2_traslocation <- 0
import$CCND2_traslocation[import$SeqFISH_T_4_14==1] <- 1
import$CCND2_traslocation[import$SeqFISH_T_14_16==1] <- 1
import$CCND2_traslocation[import$SeqFISH_T_14_20==1] <- 1

import$MMrisk_class_t_CCND2 <- ifelse(import$MMrisk_CLASS==1 & import$CCND2_traslocation ==1 , "t&1q&13+",
                                      ifelse(import$MMrisk_CLASS==1, "1q&13+_pure", import$MMrisk_allclass %>% as.character))

# import %>% select(MMrisk_class_t_CCND2, SeqFISH_T_4_14, SeqFISH_T_14_16, SeqFISH_T_14_20, MMrisk_CLASS) %>% View                                         

import$SeqFISH_T_4_14_pure <- ifelse(import$MMrisk_CLASS != 1 &  import$SeqFISH_T_4_14 ==1, 1, 0)
import$SeqFISH_T_14_16_pure <- ifelse(import$MMrisk_CLASS != 1 & import$SeqFISH_T_14_16 ==1, 1, 0)
import$SeqFISH_T_14_20_pure <- ifelse(import$MMrisk_CLASS != 1 & import$SeqFISH_T_14_20 ==1, 1, 0)


#========== Create Survival Data ============

import$OS_months <- import$ttcos / 30.5
import$PFS_months <- import$ttcpfs / 30.5

OS <- Surv(import$OS_months, import$censos)
PFS <- Surv(import$PFS_months, import$censpfs)


#=======================================================================
#========================== MULTIVARIATE ===============================
#=======================================================================


.................# ggforest2 function def #......................

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


###############################################################
################### MULTIVARIATE 03/02/22 #####################
###############################################################


#===========  MODEL 1 - multivariate MMrisk =============== 

# OS model 1
mv_R <- coxph(OS ~ OLD + SEX + 
                ISS + 
                HB_m_105 + PLT_m_150 +  
                MMrisk_allclass +
                SeqFISH_Del_1p_CDKN2C + SeqFISH_Del_TP53 + HyperDiploidy +
                IgH_translocation_type,
                 data = import )
mv_R %>% summary()

ggforest2(mv_R, main = "OS COMMpass (model 1)") 

ggsave(paste0(outpath,"FOREST_OS_COMMpass_model1.pdf"), width = 10, height = 7)


# PFS model 1
mv_R <- coxph(PFS ~ OLD + SEX + 
                ISS +
                HB_m_105 + PLT_m_150 +  
                MMrisk_allclass +
                SeqFISH_Del_1p_CDKN2C + SeqFISH_Del_TP53 + HyperDiploidy +
                IgH_translocation_type,
                data = import )
mv_R %>% summary()

ggforest2(mv_R, main = "PFS COMMpass (model 1)") 

ggsave(paste0(outpath,"FOREST_PFS_COMMpass_model1.pdf"), width = 10, height = 7)


#=========== MODEL 2 - multivariate MMrisk and CCND2 t =============== 

# OS model 2
mv_R <- coxph(OS ~ OLD + SEX + 
                ISS +
                HB_m_105 + PLT_m_150 +  
                MMrisk_class_t_CCND2 +
                SeqFISH_Del_1p_CDKN2C + SeqFISH_Del_TP53 + HyperDiploidy +
                SeqFISH_T_4_14_pure + SeqFISH_T_14_16_pure + SeqFISH_T_14_20_pure + SeqFISH_T_11_14 + SeqFISH_T_6_14,
                data = import )
mv_R %>% summary()

ggforest2(mv_R, main = "OS COMMpass (model 2)") 

ggsave(paste0(outpath,"FOREST_OS_COMMpass_model2.pdf"), width = 10, height = 7)


# PFS model 2
mv_R <- coxph(PFS ~ OLD + SEX + 
                ISS +
                HB_m_105 + PLT_m_150 +  
                MMrisk_class_t_CCND2 +
                SeqFISH_Del_1p_CDKN2C + SeqFISH_Del_TP53 + HyperDiploidy +
                SeqFISH_T_4_14_pure + SeqFISH_T_14_16_pure + SeqFISH_T_14_20_pure + SeqFISH_T_11_14 + SeqFISH_T_6_14,
                data = import )
mv_R %>% summary()

ggforest2(mv_R, main = "PFS COMMpass (model 2)") 

ggsave(paste0(outpath,"FOREST_PFS_COMMpass_model2.pdf"), width = 10, height = 7)



# #=============== ULTRA model ===================
# # OS Ultra
# mv_R <- coxph(OS ~ OLD + SEX + 
#                 HB_m_105 + PLT_m_150 +  
#                 ULTRA_Clincal + 
#                 SeqFISH_Del_1p_CDKN2C_FAF1 + SeqFISH_Del_TP53 + HyperDiploidy +
#                 IgH_translocation_type,
#                 data = import )
# mv_R %>% summary()
# 
# ggforest2(mv_R, main = "OS CoMMpass (model Ultra)") 
# 
# ggsave(paste0(outpath,"FOREST_OS_CoMMpass_ULTRA.pdf"), width = 10, height = 7)
# 
# 
# # PFS Ultra
# mv_R <- coxph(PFS ~ OLD + SEX + 
#                 HB_m_105 + PLT_m_150 +  
#                 ULTRA_Clincal +
#                 SeqFISH_Del_1p_CDKN2C_FAF1 + SeqFISH_Del_TP53 + HyperDiploidy +
#                 IgH_translocation_type,
#                 data = import )
# mv_R %>% summary()
# 
# ggforest2(mv_R, main = "PFSCoMMpass (model Ultra)") 
# 
# ggsave(paste0(outpath,"FOREST_PFS_CoMMpass_ULTRA.pdf"), width = 10, height = 7)


#=============== R_ISS model ===================
# OS Ultra
mv_R <- coxph(OS ~ OLD + SEX +
                HB_m_105 + PLT_m_150 +
                R_ISS +
                SeqFISH_Del_1p_CDKN2C + HyperDiploidy +
                SeqFISH_T_11_14,
                data = import )
mv_R %>% summary()

ggforest2(mv_R, main = "OS CoMMpass (R-ISS)")

ggsave(paste0(outpath,"FOREST_OS_CoMMpass_R-ISS.pdf"), width = 10, height = 7)


# PFS Ultra
mv_R <- coxph(PFS ~ OLD + SEX +
                HB_m_105 + PLT_m_150 +
                R_ISS +
                SeqFISH_Del_1p_CDKN2C + HyperDiploidy +
                SeqFISH_T_11_14,
                data = import )
mv_R %>% summary()

ggforest2(mv_R, main = "PFS CoMMpass (R-ISS)")

ggsave(paste0(outpath,"FOREST_PFS_CoMMpass_R-ISS.pdf"), width = 10, height = 7)

######################
# UNIVARIATE - LOOP  #
######################
#create temp functions for COX univariate
tmpfun_os <- function(x) as.formula(paste("OS",x,sep="~"))
tmpfun_pfs <- function(x) as.formula(paste("PFS",x,sep="~"))

var_co<-rev(c("OLD",
          "SEX",
          "ASCT",
          "ALBUMIN_m_3.5",
          "B2M_level",
          "LDH_level",
          "SeqFISH_Del_TP53",
          "SeqFISH_Del_1p_CDKN2C",
          "IgH_translocation_type",
          "SeqFISH_T_4_14",
          "SeqFISH_T_11_14",
          "SeqFISH_T_14_16",
          "SeqFISH_T_14_20",
          "SeqFISH_T_6_14",
          "AMP_1q_genes_all",
          
          "MMrisk_1q_all",
          "MMrisk_13_all",
          
          "MMrisk_allclass",
          
          "ISS",
          "R_ISS",
          "CCND2_traslocation",
          "MMrisk_class_t_CCND2",
          "SeqFISH_T_4_14_pure",
          "SeqFISH_T_14_16_pure",
          "SeqFISH_T_14_20_pure"))


# loop
output1<-output2<-NULL
for (i in var_co) {
  print(i)
  out1<-coxph(tmpfun_os(i), data = import) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="OS", .) 
  out2<-coxph(tmpfun_pfs(i), data = import) %>% tidy( exponentiate =T, conf.int = T) %>% cbind(surv="PFS", .)
  output1<-rbind(out1, out2)
  output2<-rbind(output1, output2)
}
output2
output<-output2[,c(2,1,3,7,8,6)]
output$p.value2<-ifelse(output$p.value<0.001, "<0.001",ifelse(output$p.value>0.05,"ns", round(output$p.value,3)))
output$sign<-ifelse(output$p.value2<0.001, "***",ifelse(output$p.value2<0.01, "**",ifelse(output$p.value2<0.05, "*", " ")))
output$p.value3<-paste(output$p.value2,output$sign, sep=" ")

output_co<-output[,c(1,2,3,4,5,9)]

names(output_co)<-c("Variable","OS/PFS", "HR","Lower CI", "Upper CI", "P-value")

sjPlot::tab_df(output_co)



#17/02 manca il commento con  l'elenco delle reference: (sotto) 

levels(import$B2M_level)[[1]]

# The reference categories are th following: 
# SEX - female, MMrisk_class_t_CCND2 - 1q&13-
# R_ISS - 2, ISS - 2, MMrisk_allclass - 1q&13-
# IgH_translocation_type - "no translocation"

# B2M_level - normal

