library(data.table)
library(tidyverse)
library(ggExtra)


batches_names <- c("HD1","HD2","SNP6")


#===================== loop over the ASCAT batches ============================

for( b in batches_names){
  
  message(b)
  
  results <- fread(paste0("data/ASCAT/batch",b, "_ASCAT.results_table.csv"))
  
  # ploidy / purity exploration and visualization
  gg <- results %>% ggplot(aes(x = ascat.output.aberrantcellfraction, y = ascat.output.ploidy, color=ascat.output.ploidy>3.5)) + 
    geom_point(alpha=0.7) + theme(legend.position = "bottom") + labs(title = "Batch 1 HD" ) 
  
  gg %>% ggExtra::ggMarginal(type = "density", groupFill = T)
  
  
  #____________ segments CN refit for WGD _______________
  
  results$WGD_event <- ifelse(results$ascat.output.ploidy>3.5, 1,0) # DEFINE a WGD event where ploidy > 3.5
  
  results$WGD_event %>% table()
  results %>% ggplot(aes(x = as.factor(WGD_event))) + geom_bar(aes(fill=as.factor(WGD_event)))
  
  # define the samples with WGD_event
  WGD_samples <- results$V1[results$WGD_event==1]
  
  
  # import segments
  segm <- fread(paste0("data/ASCAT/batch",b,"adj_fitted_and_raw_segments.tsv"))
  
  segm$name <- segm$sample %>% str_replace("_\\.CytoScanHD_Array","") 
  
  idx <- segm$sample %in% WGD_samples # define which segments belong to a WGD sample
  
 
  #______ CNvalueRaw refit _____
  # inspect
  segm %>% ggplot(aes(name, CNvalueRaw, color= as.factor(idx))) + 
    geom_boxplot(outlier.shape = NA) + ylim(c(0,6))  + theme(axis.text.x=element_blank())
  
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
    geom_boxplot(outlier.shape = NA) + ylim(c(0,6)) + theme(axis.text.x=element_blank())
  
  
  
  #____ save the refitted segments of this batch ______
  dir.create("workfiles/ASCAT_refitting")
  
  write_tsv(segm, paste0("workfiles/ASCAT_refitting/refitted_ASCAT_segments_batch",b,".tsv"))
  
  
}





