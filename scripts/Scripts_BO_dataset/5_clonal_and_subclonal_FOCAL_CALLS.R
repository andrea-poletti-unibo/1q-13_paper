
library(data.table)
library(tidyverse)
library(GenomicRanges)


# import segments from ASCAT/IGV - median cenrered (WGD correction) - factor scaled (ASCAT overestimate 1.12)
segm <- fread("results/IGV_files/IGV_Rawsegm_medianCentered_Scaled.seg") 

# visual check
segm %>% ggplot(aes(seg.mean_medAdj, weight = num.mark)) + geom_density() + xlim(-2,2) + geom_vline(xintercept = -1, colour="blue") + geom_vline(xintercept = 1, colour="red")

# import focal genes
foc_genes <- readxl::read_xlsx("input_data/focal_genes_panel_1q&13.xlsx") 

# create samples vector and SORT IT (CRITICAL!)
samples <- segm$ID %>% unique %>% sort

# create a Granges object from segments
Gsegm <- makeGRangesFromDataFrame(segm, keep.extra.columns = T)

# create a Granges object from the genes
Gfoc_genes <- makeGRangesFromDataFrame(foc_genes, keep.extra.columns = T, ignore.strand = T) 



#================ Gene-level calls algorhitm ================ 

foc_genes$focal_genes[8]
i=8

genesCallsRes <- data.frame()
gene_names <- vector()
gene_chr <- vector()

discarded_genes <- vector()

for ( i in 1:nrow(foc_genes)) {
  
  G_gene <- Gfoc_genes[i] # select a single gene
  G_gene
  wid_gene <- G_gene %>% width() # save the length of the selected gene

  overlap <- pintersect( Gsegm, G_gene) # comupte overlaps between the gene and the segment-gene positions

  res.gene <- overlap %>% as.data.frame() %>% filter(hit==TRUE) # generate the GENE-level results
  
  res.gene$percLen <- res.gene$width/wid_gene %>% round(3) # create the percenage of segment overlap with the gene
  
  res.gene <- dplyr::filter(res.gene, percLen >= 0.10) # exclude those segments that overlap the gene for less then 10%
  
  res.gene <- res.gene %>% arrange(ID, desc(seg.mean_medAdj)) # ! SUPER IMPORTANT ORDERING: per sample and then per descending CN change -> the first entry per sample will be the most changed in CN!
  
  res.gene <- dplyr::distinct(res.gene, ID, .keep_all = TRUE ) # KEEP only the first entry for each sample (this will end up always with a number of rows equal to number of samples) 
  
  CNvalues<- res.gene$seg.mean_medAdj # extract the CN values 
  
  # SAVE data ONLY if the calls are complete: there should be a CN call for EACH sample
  if (length(CNvalues) == length(samples) ) { 
    
    genesCallsRes <- rbind( genesCallsRes, CNvalues) # save the CN values
    gene_names <- append(gene_names, G_gene$focal_genes) # save the gene name
    gene_chr <- append(gene_chr, as.character(G_gene@seqnames)) # save the gene chr
    
  } else { 
    discarded_genes <- append(discarded_genes, G_gene$focal_genes )} # save the discarded incomplete gene
  
}

names(genesCallsRes) <- samples
rownames(genesCallsRes) <- gene_names


genesCallsFinal <- t(genesCallsRes) # traspose the results



#================ compute the 1/0 calls ================
library(gplots)


dir.create("workfiles/FOCAL_CALLS")

#______ call MAJOR (clonality > 0.5) _______

# AMPS
CALLS_AMP_MAJOR <- apply(as.matrix(genesCallsFinal), c(1,2), function(x) if(x>0.5) 1 else 0) %>% as.data.frame()
names(CALLS_AMP_MAJOR) <- paste0("AMP_maj_focal_", names(CALLS_AMP_MAJOR))

heatmap.2(as.matrix(CALLS_AMP_MAJOR), col=c("grey","blue"), trace= "none")

CALLS_AMP_MAJOR_sel <- CALLS_AMP_MAJOR %>% select(`AMP_maj_focal_MCL1`:`AMP_maj_focal_MYC`)

write.table(CALLS_AMP_MAJOR_sel,"workfiles/FOCAL_CALLS/AMP_CALLS_major_focal.txt", sep = "\t", quote = F, row.names = T)

# DELS
CALLS_DEL_MAJOR <- apply(as.matrix(genesCallsFinal), c(1,2), function(x) if(x< -0.5) 1 else 0) %>% as.data.frame()
names(CALLS_DEL_MAJOR) <- paste0("DEL_maj_focal_", names(CALLS_DEL_MAJOR))

heatmap.2(as.matrix(CALLS_DEL_MAJOR), col=c("grey","red"), trace= "none")

CALLS_DEL_MAJOR_sel <- CALLS_DEL_MAJOR %>% select(`DEL_maj_focal_FAM46C`:`DEL_maj_focal_TP53`)
write.table(CALLS_DEL_MAJOR_sel,"workfiles/FOCAL_CALLS/DEL_CALLS_major_focal.txt", sep = "\t", quote = F, row.names = T)


#______ call CLONAL (clonality > 0.9) _______

# AMPS
CALLS_AMP_CLONAL <- apply(as.matrix(genesCallsFinal), c(1,2), function(x) if(x>0.9) 1 else 0) %>% as.data.frame()
names(CALLS_AMP_CLONAL) <- paste0("AMP_clon_focal_", names(CALLS_AMP_CLONAL))

heatmap.2(as.matrix(CALLS_AMP_CLONAL), col=c("grey","blue"), trace= "none")

CALLS_AMP_CLONAL_sel <- CALLS_AMP_CLONAL %>% select(`AMP_clon_focal_MCL1`:`AMP_clon_focal_MYC`)
write.table(CALLS_AMP_CLONAL_sel, "workfiles/FOCAL_CALLS/AMP_CALLS_clonal_focal.txt", sep = "\t", quote = F, row.names = T)

# DELS
CALLS_DEL_CLONAL <- apply(as.matrix(genesCallsFinal), c(1,2), function(x) if(x< -0.9) 1 else 0) %>% as.data.frame()
names(CALLS_DEL_CLONAL) <- paste0("DEL_clon_focal_", names(CALLS_DEL_CLONAL))

heatmap.2(as.matrix(CALLS_DEL_CLONAL), col=c("grey","red"), trace= "none")

CALLS_DEL_CLONAL_sel <- CALLS_DEL_CLONAL %>% select(`DEL_clon_focal_FAM46C`:`DEL_clon_focal_TP53`)
write.table(CALLS_DEL_CLONAL_sel,"workfiles/FOCAL_CALLS/DEL_CALLS_clonal_focal.txt", sep = "\t", quote = F, row.names = T)



#______ call ALL (clonality > 0.1) _______

# AMPS
CALLS_AMP_ALL <- apply(as.matrix(genesCallsFinal), c(1,2), function(x) if(x>0.1) 1 else 0) %>% as.data.frame()
names(CALLS_AMP_ALL) <- paste0("AMP_all_focal_", names(CALLS_AMP_ALL))

heatmap.2(as.matrix(CALLS_AMP_ALL), col=c("grey","blue"), trace= "none")

CALLS_AMP_ALL_sel <- CALLS_AMP_ALL %>% select(`AMP_all_focal_MCL1`:`AMP_all_focal_MYC`)
write.table(CALLS_AMP_ALL_sel, "workfiles/FOCAL_CALLS/AMP_CALLS_all_focal.txt", sep = "\t", quote = F, row.names = T)

# DELS
CALLS_DEL_ALL <- apply(as.matrix(genesCallsFinal), c(1,2), function(x) if(x< -0.1) 1 else 0) %>% as.data.frame()
names(CALLS_DEL_ALL) <- paste0("DEL_all_focal_", names(CALLS_DEL_ALL))

heatmap.2(as.matrix(CALLS_DEL_ALL), col=c("grey","red"), trace= "none")

CALLS_DEL_ALL_sel <- CALLS_DEL_ALL %>% select(`DEL_all_focal_FAM46C`:`DEL_all_focal_TP53`)
write.table(CALLS_DEL_ALL_sel,"workfiles/FOCAL_CALLS/DEL_CALLS_all_focal.txt", sep = "\t", quote = F, row.names = T)



