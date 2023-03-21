library(data.table)
library(tidyverse)
library(ggExtra)


################## BATCH 1 HD ####################

setwd("D:/analisi_in_corso/risk_1q_13/ASCAT/batchHD1/")

results <- fread("batchHD1_ASCAT.results_table.csv")

(results %>% ggplot(aes(x = ascat.output.aberrantcellfraction, y = ascat.output.ploidy, color=ascat.output.ploidy>3.5)) + 
  geom_point(alpha=0.7) + theme(legend.position = "bottom") + labs(title = "Batch 1 HD" ) ) %>% ggExtra::ggMarginal(type = "density", groupFill = T)


results %>% ggplot(aes(x = ascat.output.aberrantcellfraction, y = ascat.output.ploidy, label= V1)) + 
  geom_point() + geom_text(size=2, nudge_y = 0.04)

#____________ segments CN refit for WGD _______________

results$WGD_event <- ifelse(results$ascat.output.ploidy>3.5, 1,0) # DEFINE a WGD event where ploidy > 3.5
results$WGD_event %>% table()
21/190

results %>% ggplot(aes(x = as.factor(WGD_event))) + geom_bar(aes(fill=as.factor(WGD_event)))

WGD_samples <- results$V1[results$WGD_event==1]

segm <- fread("batchHD1adj_fitted_and_raw_segments.tsv")

segm$name <- segm$sample %>% str_replace("_\\.CytoScanHD_Array","") 

idx <- segm$sample %in% WGD_samples # define which segments belong to a WGD sample

#_____ CNvalue refit_____
# inspect
segm %>% ggplot(aes(name, CNvalue, color= as.factor(idx))) + 
  geom_boxplot(outlier.shape = NA) + ylim(c(0,6)) + theme(axis.text.x = element_text(angle = 90))

# refitting
CNvalue_refit <- vector() 
for(i in 1:length(idx)){
  if (idx[i]==TRUE) {
    CNvalue_refit <- append(CNvalue_refit, segm$CNvalue[i]/2)
  } else {
    CNvalue_refit <- append(CNvalue_refit, segm$CNvalue[i])
  }
}
segm$CNvalue_refit <- CNvalue_refit

# Check the effect of refitting
segm %>% ggplot(aes(name, CNvalue_refit, color= as.factor(idx))) + 
  geom_boxplot(outlier.shape = NA) + ylim(c(0,6)) + theme(axis.text.x = element_text(angle = 90))


#______ CNvalueRaw refit _____
# inspect
segm %>% ggplot(aes(name, CNvalueRaw, color= as.factor(idx))) + 
  geom_boxplot(outlier.shape = NA) + ylim(c(0,6)) + theme(axis.text.x = element_text(angle = 90))

# refitting
CNvalueRaw_refit <- vector() 
for(i in 1:length(idx)){
  if (idx[i]==TRUE) {
    CNvalueRaw_refit <- append(CNvalueRaw_refit, segm$CNvalueRaw[i]/2)
  } else {
    CNvalueRaw_refit <- append(CNvalueRaw_refit, segm$CNvalueRaw[i])
  }
}
segm$CNvalueRaw_refit <- CNvalueRaw_refit

# Check the effect of refitting
segm %>% ggplot(aes(name, CNvalueRaw_refit, color= as.factor(idx))) + 
  geom_boxplot(outlier.shape = NA) + ylim(c(0,6)) + theme(axis.text.x = element_text(angle = 90))



#____ save the refitted segments of this batch ______
dir.create("D:/analisi_in_corso/risk_1q_13/ASCAT/Refitted_ASCAT_segm")
write_tsv(segm, "D:/analisi_in_corso/risk_1q_13/ASCAT/Refitted_ASCAT_segm/refitted_ASCAT_segments_batchHD1.tsv")


################## BATCH 2 HD ####################

setwd("D:/analisi_in_corso/risk_1q_13/ASCAT/batchHD2/")

results <- fread("batchHD2_ASCAT.results_table.csv")

(results %>% ggplot(aes(x = ascat.output.aberrantcellfraction, y = ascat.output.ploidy, color=ascat.output.ploidy>3.5)) + 
    geom_point(alpha=0.7) + theme(legend.position = "bottom")+ labs(title = "Batch 2 HD" )) %>% ggExtra::ggMarginal(type = "density", groupFill = T)


results %>% ggplot(aes(x = ascat.output.aberrantcellfraction, y = ascat.output.ploidy, label= V1)) + 
  geom_point() + geom_text(size=2, nudge_y = 0.04)

#____________ segments CN refit for WGD _______________

results$WGD_event <- ifelse(results$ascat.output.ploidy>3.5, 1,0) # DEFINE a WGD event where ploidy > 3.5
results$WGD_event %>% table()
31/194

results %>% ggplot(aes(x = as.factor(WGD_event))) + geom_bar(aes(fill=as.factor(WGD_event)))

WGD_samples <- results$V1[results$WGD_event==1]

segm <- fread("batchHD2adj_fitted_and_raw_segments.tsv")

segm$name <- segm$sample %>% str_replace("_\\.CytoScanHD_Array","") 

idx <- segm$sample %in% WGD_samples # define which segments belong to a WGD sample

#_____ CNvalue refit_____
# inspect
segm %>% ggplot(aes(name, CNvalue, color= as.factor(idx))) + 
  geom_boxplot(outlier.shape = NA) + ylim(c(0,6)) + theme(axis.text.x = element_text(angle = 90))

# refitting
CNvalue_refit <- vector() 
for(i in 1:length(idx)){
  if (idx[i]==TRUE) {
    CNvalue_refit <- append(CNvalue_refit, segm$CNvalue[i]/2)
  } else {
    CNvalue_refit <- append(CNvalue_refit, segm$CNvalue[i])
  }
}
segm$CNvalue_refit <- CNvalue_refit

# Check the effect of refitting
segm %>% ggplot(aes(name, CNvalue_refit, color= as.factor(idx))) + 
  geom_boxplot(outlier.shape = NA) + ylim(c(0,6)) + theme(axis.text.x = element_text(angle = 90))


#______ CNvalueRaw refit _____
# inspect
segm %>% ggplot(aes(name, CNvalueRaw, color= as.factor(idx))) + 
  geom_boxplot(outlier.shape = NA) + ylim(c(0,6)) + theme(axis.text.x = element_text(angle = 90))

# refitting
CNvalueRaw_refit <- vector() 
for(i in 1:length(idx)){
  if (idx[i]==TRUE) {
    CNvalueRaw_refit <- append(CNvalueRaw_refit, segm$CNvalueRaw[i]/2)
  } else {
    CNvalueRaw_refit <- append(CNvalueRaw_refit, segm$CNvalueRaw[i])
  }
}
segm$CNvalueRaw_refit <- CNvalueRaw_refit

# Check the effect of refitting
segm %>% ggplot(aes(name, CNvalueRaw_refit, color= as.factor(idx))) + 
  geom_boxplot(outlier.shape = NA) + ylim(c(0,6)) + theme(axis.text.x = element_text(angle = 90))

#____ save the refitted segments of this batch ______
write_tsv(segm, "D:/analisi_in_corso/risk_1q_13/ASCAT/Refitted_ASCAT_segm/refitted_ASCAT_segments_batchHD2.tsv")



################## BATCH SNP 6 ####################

setwd("D:/analisi_in_corso/risk_1q_13/ASCAT/batchSNP6/")

results <- fread("batchSNP6_ASCAT.results_table.csv")

(results %>% ggplot(aes(x = ascat.output.aberrantcellfraction, y = ascat.output.ploidy, color=ascat.output.ploidy>3.5)) + 
    geom_point(alpha=0.7) + theme(legend.position = "bottom") + labs(title = "Batch SNP6" ) ) %>% ggExtra::ggMarginal(type = "density", groupFill = T)


results %>% ggplot(aes(x = ascat.output.aberrantcellfraction, y = ascat.output.ploidy, label= V1)) + 
  geom_point() + geom_text(size=2, nudge_y = 0.04)

#____________ segments CN refit for WGD _______________

results$WGD_event <- ifelse(results$ascat.output.ploidy>3.5, 1,0) # DEFINE a WGD event where ploidy > 3.5
results$WGD_event %>% table()
23/131

results %>% ggplot(aes(x = as.factor(WGD_event))) + geom_bar(aes(fill=as.factor(WGD_event)))

WGD_samples <- results$V1[results$WGD_event==1]

segm <- fread("batchSNP6adj_fitted_and_raw_segments.tsv")

segm$name <- segm$sample %>% str_replace("_\\.GenomeWideSNP_6\\.","") 

idx <- segm$sample %in% WGD_samples # define which segments belong to a WGD sample

#_____ CNvalue refit_____
# inspect
segm %>% ggplot(aes(name, CNvalue, color= as.factor(idx))) + 
  geom_boxplot(outlier.shape = NA) + ylim(c(0,6)) + theme(axis.text.x = element_text(angle = 90))

# refitting
CNvalue_refit <- vector() 
for(i in 1:length(idx)){
  if (idx[i]==TRUE) {
    CNvalue_refit <- append(CNvalue_refit, segm$CNvalue[i]/2)
  } else {
    CNvalue_refit <- append(CNvalue_refit, segm$CNvalue[i])
  }
}
segm$CNvalue_refit <- CNvalue_refit

# Check the effect of refitting
segm %>% ggplot(aes(name, CNvalue_refit, color= as.factor(idx))) + 
  geom_boxplot(outlier.shape = NA) + ylim(c(0,6)) + theme(axis.text.x = element_text(angle = 90))


#______ CNvalueRaw refit _____
# inspect
segm %>% ggplot(aes(name, CNvalueRaw, color= as.factor(idx))) + 
  geom_boxplot(outlier.shape = NA) + ylim(c(0,6)) + theme(axis.text.x = element_text(angle = 90))

# refitting
CNvalueRaw_refit <- vector() 
for(i in 1:length(idx)){
  if (idx[i]==TRUE) {
    CNvalueRaw_refit <- append(CNvalueRaw_refit, segm$CNvalueRaw[i]/2)
  } else {
    CNvalueRaw_refit <- append(CNvalueRaw_refit, segm$CNvalueRaw[i])
  }
}
segm$CNvalueRaw_refit <- CNvalueRaw_refit

# Check the effect of refitting
segm %>% ggplot(aes(name, CNvalueRaw_refit, color= as.factor(idx))) + 
  geom_boxplot(outlier.shape = NA) + ylim(c(0,6)) + theme(axis.text.x = element_text(angle = 90))

#____ save the refitted segments of this batch ______
write_tsv(segm, "D:/analisi_in_corso/risk_1q_13/ASCAT/Refitted_ASCAT_segm/refitted_ASCAT_segments_batchSNP6.tsv")


# deprecated ################# BATCH Failed HD #########################

# setwd("D:/analisi_in_corso/risk_1q_13/ASCAT/batchHDF/")
# 
# results <- fread("batchHDF_ASCAT.results_table.csv")
# 
# (results %>% ggplot(aes(x = ascat.output.aberrantcellfraction, y = ascat.output.ploidy, color=ascat.output.ploidy>3.5)) + 
#     geom_point(alpha=0.7) + theme(legend.position = "bottom")) %>% ggExtra::ggMarginal(type = "density", groupFill = T)
# 
# 
# results %>% ggplot(aes(x = ascat.output.aberrantcellfraction, y = ascat.output.ploidy, label= V1)) + 
#   geom_point() + geom_text(size=2, nudge_y = 0.04)
# 
# #____________ segments CN refit for WGD _______________
# 
# results$WGD_event <- ifelse(results$ascat.output.ploidy>3.5, 1,0) # DEFINE a WGD event where ploidy > 3.5
# results$WGD_event %>% table()
# 8/36
# 
# results %>% ggplot(aes(x = as.factor(WGD_event))) + geom_bar(aes(fill=as.factor(WGD_event)))
# 
# WGD_samples <- results$V1[results$WGD_event==1]
# 
# segm <- fread("batchHDFadj_fitted_and_raw_segments.tsv")
# 
# segm$name <- segm$sample %>% str_replace("_\\.CytoScanHD_Array","") 
# 
# idx <- segm$sample %in% WGD_samples # define which segments belong to a WGD sample
# 
# #_____ CNvalue refit_____
# # inspect
# segm %>% ggplot(aes(name, CNvalue, color= as.factor(idx))) + 
#   geom_boxplot(outlier.shape = NA) + ylim(c(0,6)) + theme(axis.text.x = element_text(angle = 90))
# 
# # refitting
# CNvalue_refit <- vector() 
# for(i in 1:length(idx)){
#   if (idx[i]==TRUE) {
#     CNvalue_refit <- append(CNvalue_refit, segm$CNvalue[i]/2)
#   } else {
#     CNvalue_refit <- append(CNvalue_refit, segm$CNvalue[i])
#   }
# }
# segm$CNvalue_refit <- CNvalue_refit
# 
# # Check the effect of refitting
# segm %>% ggplot(aes(name, CNvalue_refit, color= as.factor(idx))) + 
#   geom_boxplot(outlier.shape = NA) + ylim(c(0,6)) + theme(axis.text.x = element_text(angle = 90))
# 
# 
# #______ CNvalueRaw refit _____
# # inspect
# segm %>% ggplot(aes(name, CNvalueRaw, color= as.factor(idx))) + 
#   geom_boxplot(outlier.shape = NA) + ylim(c(0,6)) + theme(axis.text.x = element_text(angle = 90))
# 
# # refitting
# CNvalueRaw_refit <- vector() 
# for(i in 1:length(idx)){
#   if (idx[i]==TRUE) {
#     CNvalueRaw_refit <- append(CNvalueRaw_refit, segm$CNvalueRaw[i]/2)
#   } else {
#     CNvalueRaw_refit <- append(CNvalueRaw_refit, segm$CNvalueRaw[i])
#   }
# }
# segm$CNvalueRaw_refit <- CNvalueRaw_refit
# 
# # Check the effect of refitting
# segm %>% ggplot(aes(name, CNvalueRaw_refit, color= as.factor(idx))) + 
#   geom_boxplot(outlier.shape = NA) + ylim(c(0,6)) + theme(axis.text.x = element_text(angle = 90))
# 
# #____ save the refitted segments of this batch ______
# write_tsv(segm, "D:/analisi_in_corso/risk_1q_13/ASCAT/Refitted_ASCAT_segm/refitted_ASCAT_segments_batchHDF.tsv")
# 
# 
# 
