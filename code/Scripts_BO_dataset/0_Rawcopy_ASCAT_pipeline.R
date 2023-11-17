# library(devtools)

# # install_bitbucket('n0s3n/rawcopy/rawcopyPackage')
# # install_github('VanLoo-lab/ascat/ASCAT')

# library(rawcopy)
# library(ASCAT)
# library(data.table)
# library(tidyverse)
# library(readr)


# RAWCOPY_INPUT <- "set/your/path/here"
# RAWCOPY_OUTPUT <- "set/your/path/here"

# ASCAT_dir <- "set/your/path/here"


# ################################# RAWCOPY #################################

# # 1) run rawcopy for every CEL file in INPUT directory

# rawcopy(CELfiles.or.directory = RAWCOPY_INPUT, # INPUT: directory with CEL files
#         referenceFile = NULL,
#         wavinessFile = NULL, 
#         outdir = RAWCOPY_OUTPUT, # OUTPUT: directory with rawcopy results folders
#         cores = 8, 
#         segmentation = TRUE, # Whether to produce segment table files for each sample. Currently single-threaded.
#         writeData = FALSE, # Whether to write log-ratio and B-allele ratio to text files.
#         writeSummary = TRUE, # Whether to summarize sample data in a table file.
#         allSNPs = FALSE, # Whether to report B-allele ratio for SNPs of poor quality (based on reference material only).
#         resume = FALSE, 
#         segment.type = "PSCBS", 
#         segment.alpha = 10^-7, # Sets the segmentation significance threshold. A lower value will result in fewer segments.
#         seed = NULL, 
#         MDS = TRUE)


# # 2)  data preparation for ASCAT analysis

# CELfiles <- list.files(RAWCOPY_INPUT, pattern = ".CEL", full.names = F) %>% str_remove(".CEL")

# idx<- grep("CytoScanHD", CELfiles) # index for HD
# filesHD<-CELfiles[idx]

# idx<- grep("SNP_6", CELfiles) # index for 6.0
# filesSNP6<-CELfiles[idx]


# files1<-filesHD[1:190] # BATCH 1
# files2<-filesHD[190:383] # BATCH 2
# files3<-filesSNP6 # BATCH 3

# batches <- list(batch_1_HD = files1, 
#                 batch_2_HD = files2, 
#                 batch_3_SNP6 = files3)


# # loop for every batch
# for(j in 1:length(batches)){

#   bc <- batches[[j]]
#   bc.name <- names(batches[j])
  
#   # set base working directory (Rawcopy folders dir)
#   basewd<- RAWCOPY_OUTPUT
  
#   # LOOP to extract logR and BAF from the .Rdata rawcopy output of each sample
#   result.logR<-list()
#   result.BAF<-list()
  
#   for (i in bc) { 
#     t1<-Sys.time()
    
#     print(i)
#     setwd(paste0(basewd,"/",i)) # change wd to a new sample
#     load("rawcopy.Rdata") #load the sample specific .Rdata file
#     result.logR[[i]]<-probes.txt$Value #add the sample logR to the total list
#     result.BAF[[i]]<-probes.txt$B.Allele.Frequency # add the sample BAF to the total list
    
#     t2<-Sys.time()
#     print(t2-t1)
#     }

#   # transforming the list to a data frame
#   final.logR<-as.data.frame(result.logR) 
#   final.BAF<-as.data.frame(result.BAF)
  
#   # add columns with probe info
#   complete.logR<-base::cbind(probes.txt[,1:3],final.logR)
#   rm(final.logR)
  
#   complete.BAF<- base::cbind(probes.txt[,1:3],final.BAF)
#   rm(final.BAF)
  
#   # removing "chr" in Chromosome column
#   complete.logR$Chromosome<- gsub("chr","", complete.logR$Chromosome)
#   complete.BAF$Chromosome<- gsub("chr","", complete.BAF$Chromosome)
  
#   # renaming columns in ASCAT format
#   colnames(complete.logR)[1:3]<-c("probe","chrs","pos")
#   colnames(complete.BAF)[1:3]<-c("probe","chrs","pos")
  
#   # remove sex chromosomes and Mitochondrial 
#   LOGR_noXYM <- filter(complete.logR, chrs != "M",chrs != "X", chrs !="Y")
#   BAF_noXYM <- filter(complete.BAF, chrs != "M",chrs != "X", chrs !="Y")
  
  
#   # exporting the files
  
#   dir.create(paste0(ASCAT_dir, bc.name))
  
#   setwd(paste0(ASCAT_dir, bc.name))
  
#   write_tsv( LOGR_noXYM, paste0("ascatReady.noXYM_logR_", project, "_", bc.name,".txt"))
#   write_tsv( BAF_noXYM, paste0("ascatReady.noXYM_BAF_", project, "_", bc.name,".txt"))

# }


# ################################# ASCAT #################################

# library(foreach)
# library(doParallel)
# library(parallel)


# CORES = detectCores() - 2
# registerDoParallel(CORES) 


# path <- ASCAT_dir

# batches <- list.dirs(path, recursive = F, full.names = F)
  
# foreach(i=batches, .packages = c("ASCAT", "data.table", "tidyverse", "readr") ) %dopar% {
  
#   print(i)

#   type <- ifelse(str_detect(i, "HD"), "AffyCytoScanHD", "AffySNP6")
  
#   setwd(paste0(path, i))
  
  
#   #_________ PART 1: set up dir tree and files _________
  
#   #step 1: identify multi-sample logR and BAF files (pre-generated and already present in the working directory) 
#   file_logR <-list.files(pattern="LogR")
#   file_BAF  <-list.files(pattern="BAF")
  
  
#   #step 2: creation of sub-directories for the different GISTIC phases outputs
#   dir.create("1.rawdata_plots")
#   dir.create("2.Sep_plots")
#   dir.create("3.segment_data")
#   dir.create("4.segmented_plots")
#   dir.create("5.ASCAT_output_plots")
  
  
#   #_________ PART 2: running ASCAT_________
#   #step 0: loading data
#   ascat.bc = ascat.loadData(file_logR, 
#                             file_BAF, 
#                             chrs = c(paste0("chr",1:22),"chrX","chrY","chrM"))
  
  
#   #step 1: plotting raw logR and BAF data 
#   ascat.plotRawData(ascat.bc,
#                     img.dir = "1.rawdata_plots")
  
#   #step 2: predict genotypes areas from BAF data
#   gg<-ascat.predictGermlineGenotypes(ascat.bc, 
#                                      platform = type, # AffyCytoScanHD or AffySNP6
#                                      img.dir = "2.Sep_plots")
  
#   #step 3: running ASPCF segmentation on both logR and BAF
#   ascat.bc = ascat.aspcf(ascat.bc, 
#                          ascat.gg=gg, 
#                          penalty=70, #default = 25 ( bigger => less segments)
#                          out.dir= "3.segment_data") 
  
#   #step 4: plotting segmented data
#   ascat.plotSegmentedData(ascat.bc,
#                           img.dir="4.segmented_plots")
  
#   #step 5: running ASCAT analysis and plotting outputs
#   ascat.output = ascat.runAscat(ascat.bc,
#                                 gamma=0.55, # default= 0.55
#                                 y_limit=5, # default= 5
#                                 circos = NA,
#                                 rho_manual = NA, # rho= aberrant cell fraction
#                                 psi_manual = NA, # psi= ploidy
#                                 img.dir = "5.ASCAT_output_plots")
  
  
#   #_________ part 3: exporting data results_________
  
#   # step 1: export a csv table of results
#   exportRes<- data.frame(ascat.output$goodnessOfFit,
#                          ascat.output$aberrantcellfraction,
#                          ascat.output$ploidy,
#                          ascat.output$psi)
  
#   write.csv(exportRes, paste0( i, "_ASCAT.results_table.csv"))
  
#   #step 2: export the lists of failed arrays and non-aberrant arrays
#   write.csv(ascat.output$failedarrays, paste0( i,"_failed_arrays_list.csv"))
#   write.csv(ascat.output$nonaberrantarrays, paste0( i,"_nonAberrant_arrays_list.csv"))
  
#   #step 3: export a tsv file with the segments for every sample (adjusted-fitted and adjusted-raw)
#   segmentsAdj<-ascat.output$segments
#   segmentsAdj$CNvalue<- segmentsAdj$nMajor+segmentsAdj$nMinor
#   segmentsAdj$length<- segmentsAdj$endpos - segmentsAdj$startpos
#   write_tsv(segmentsAdj, paste0( i,"_adj_fitted_segments.tsv"))
  
#   segmentsAdjRaw<-ascat.output$segments_raw
#   segmentsAdjRaw$CNvalue<- segmentsAdjRaw$nMajor + segmentsAdjRaw$nMinor
#   segmentsAdjRaw$CNvalueRaw<- segmentsAdjRaw$nAraw + segmentsAdjRaw$nBraw
#   segmentsAdjRaw$length<- segmentsAdjRaw$endpos - segmentsAdjRaw$startpos
#   write_tsv(segmentsAdjRaw, paste0( i,"adj_fitted_and_raw_segments.tsv"))
  
#   # export ascat.output as a .Rdata single object - use readRDS() function to restore/load it in a new R session
#   saveRDS(ascat.output, paste0( i,"ascat_output.RData"))

# }

# stopImplicitCluster()
