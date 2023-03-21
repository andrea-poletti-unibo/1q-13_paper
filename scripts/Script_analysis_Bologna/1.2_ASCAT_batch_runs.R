
########################## BATCH HD1 ##############################
library(ASCAT)


directories <- c("batchHD1")

i = directories[1]

#_________PART 1: set up dir tree and files_________

#step 0: choose working directory
# setwd(choose.dir()) 
setwd(paste0("D:/analisi_in_corso/risk_1q_13/ASCAT/", i))

#step 1: identify multi-sample logR and BAF files (pre-generated and already present in the working directory) 
file_logR<-list.files(pattern="logR")
file_BAF <-list.files(pattern="BAF")

#step 2: creation of sub-directories for the different GISTIC phases outputs
dir.create("1.rawdata_plots")
dir.create("2.Sep_plots")
dir.create("3.segment_data")
dir.create("4.segmented_plots")
dir.create("5.ASCAT_output_plots")


#_________PART 2: running ASCAT_________


#step 0: loading data
ascat.bc = ascat.loadData(file_logR, file_BAF, 
                          chrs = c(1:22))


#step 1: plotting raw logR and BAF data 
ascat.plotRawData(ascat.bc,
                  img.dir = "1.rawdata_plots")


#step 2: predict genotypes areas from BAF data
gg<-ascat.predictGermlineGenotypes(ascat.bc, 
                                   platform = "AffySNP6",
                                   img.dir = "2.Sep_plots")


#step 3: running ASPCF segmentation on both logR and BAF
ascat.bc = ascat.aspcf(ascat.bc, 
                       ascat.gg=gg, 
                       penalty=70, #default = 25 ( bigger => less segments)
                       out.dir= "3.segment_data") 


#step 4: plotting segmented data
ascat.plotSegmentedData(ascat.bc,
                        img.dir="4.segmented_plots")


#step 5: running ASCAT analysis and plotting outputs
ascat.output = ascat.runAscat(ascat.bc,
                              gamma=0.55, # default= 0.55
                              y_limit=5, # default= 5
                              circos = NA,
                              rho_manual = NA, # rho= aberrant cell fraction
                              psi_manual = NA, # psi= ploidy
                              img.dir = "5.ASCAT_output_plots")


#_________ part 3: exporting data results_________

# #visualize results
# ascat.output$aberrantcellfraction
# ascat.output$ploidy
# ascat.output$psi
# ascat.output$goodnessOfFit
# ascat.output$failedarrays
# ascat.output$nonaberrantarrays

# step 1: export a csv table of results
exportRes<- data.frame(ascat.output$goodnessOfFit,
                       ascat.output$aberrantcellfraction,
                       ascat.output$ploidy,
                       ascat.output$psi)

write.csv(exportRes, paste0( i, "_ASCAT.results_table.csv"))

#step 2: export the lists of failed arrays and non-aberrant arrays

write.csv(ascat.output$failedarrays, paste0( i,"_failed_arrays_list.csv"))
write.csv(ascat.output$nonaberrantarrays, paste0( i,"_nonAberrant_arrays_list.csv"))

#step 3: export a tsv file with the segments for every sample (adjusted-fitted and adjusted-raw)
library(readr)

segmentsAdj<-ascat.output$segments
segmentsAdj$CNvalue<- segmentsAdj$nMajor+segmentsAdj$nMinor
segmentsAdj$length<- segmentsAdj$endpos - segmentsAdj$startpos
write_tsv(segmentsAdj, paste0( i,"_adj_fitted_segments.tsv"))

segmentsAdjRaw<-ascat.output$segments_raw
segmentsAdjRaw$CNvalue<- segmentsAdjRaw$nMajor + segmentsAdjRaw$nMinor
segmentsAdjRaw$CNvalueRaw<- segmentsAdjRaw$nAraw + segmentsAdjRaw$nBraw
segmentsAdjRaw$length<- segmentsAdjRaw$endpos - segmentsAdjRaw$startpos
write_tsv(segmentsAdjRaw, paste0( i,"adj_fitted_and_raw_segments.tsv"))

# export ascat.output as a .Rdata single object - use readRDS() function to restore/load it in a new R session
saveRDS(ascat.output, paste0( i,"ascat_output.RData"))

# my_output<-readRDS("ascat_output.RData")



########################## BATCH HD2 ##############################
library(ASCAT)


directories <- c("batchHD2")

i = directories[1]

#_________PART 1: set up dir tree and files_________

#step 0: choose working directory
# setwd(choose.dir()) 
setwd(paste0("D:/analisi_in_corso/risk_1q_13/ASCAT/", i))

#step 1: identify multi-sample logR and BAF files (pre-generated and already present in the working directory) 
file_logR<-list.files(pattern="logR")
file_BAF <-list.files(pattern="BAF")

#step 2: creation of sub-directories for the different GISTIC phases outputs
dir.create("1.rawdata_plots")
dir.create("2.Sep_plots")
dir.create("3.segment_data")
dir.create("4.segmented_plots")
dir.create("5.ASCAT_output_plots")


#_________PART 2: running ASCAT_________


#step 0: loading data
ascat.bc = ascat.loadData(file_logR, file_BAF, 
                          chrs = c(1:22))


#step 1: plotting raw logR and BAF data 
ascat.plotRawData(ascat.bc,
                  img.dir = "1.rawdata_plots")


#step 2: predict genotypes areas from BAF data
gg<-ascat.predictGermlineGenotypes(ascat.bc, 
                                   platform = "AffySNP6",
                                   img.dir = "2.Sep_plots")


#step 3: running ASPCF segmentation on both logR and BAF
ascat.bc = ascat.aspcf(ascat.bc, 
                       ascat.gg=gg, 
                       penalty=70, #default = 25 ( bigger => less segments)
                       out.dir= "3.segment_data") 


#step 4: plotting segmented data
ascat.plotSegmentedData(ascat.bc,
                        img.dir="4.segmented_plots")


#step 5: running ASCAT analysis and plotting outputs
ascat.output = ascat.runAscat(ascat.bc,
                              gamma=0.55, # default= 0.55
                              y_limit=5, # default= 5
                              circos = NA,
                              rho_manual = NA, # rho= aberrant cell fraction
                              psi_manual = NA, # psi= ploidy
                              img.dir = "5.ASCAT_output_plots")


#_________ part 3: exporting data results_________

# #visualize results
# ascat.output$aberrantcellfraction
# ascat.output$ploidy
# ascat.output$psi
# ascat.output$goodnessOfFit
# ascat.output$failedarrays
# ascat.output$nonaberrantarrays

# step 1: export a csv table of results
exportRes<- data.frame(ascat.output$goodnessOfFit,
                       ascat.output$aberrantcellfraction,
                       ascat.output$ploidy,
                       ascat.output$psi)

write.csv(exportRes, paste0( i, "_ASCAT.results_table.csv"))

#step 2: export the lists of failed arrays and non-aberrant arrays

write.csv(ascat.output$failedarrays, paste0( i,"_failed_arrays_list.csv"))
write.csv(ascat.output$nonaberrantarrays, paste0( i,"_nonAberrant_arrays_list.csv"))

#step 3: export a tsv file with the segments for every sample (adjusted-fitted and adjusted-raw)
library(readr)

segmentsAdj<-ascat.output$segments
segmentsAdj$CNvalue<- segmentsAdj$nMajor+segmentsAdj$nMinor
segmentsAdj$length<- segmentsAdj$endpos - segmentsAdj$startpos
write_tsv(segmentsAdj, paste0( i,"_adj_fitted_segments.tsv"))

segmentsAdjRaw<-ascat.output$segments_raw
segmentsAdjRaw$CNvalue<- segmentsAdjRaw$nMajor + segmentsAdjRaw$nMinor
segmentsAdjRaw$CNvalueRaw<- segmentsAdjRaw$nAraw + segmentsAdjRaw$nBraw
segmentsAdjRaw$length<- segmentsAdjRaw$endpos - segmentsAdjRaw$startpos
write_tsv(segmentsAdjRaw, paste0( i,"adj_fitted_and_raw_segments.tsv"))

# export ascat.output as a .Rdata single object - use readRDS() function to restore/load it in a new R session
saveRDS(ascat.output, paste0( i,"ascat_output.RData"))

# my_output<-readRDS("ascat_output.RData")



########################## BATCH SNP6 #######################
library(ASCAT)


directories <- c("batchSNP6")

i = directories[1]

#_________PART 1: set up dir tree and files_________

#step 0: choose working directory
# setwd(choose.dir()) 
setwd(paste0("D:/analisi_in_corso/risk_1q_13/ASCAT/", i))

#step 1: identify multi-sample logR and BAF files (pre-generated and already present in the working directory) 
file_logR<-list.files(pattern="logR")
file_BAF <-list.files(pattern="BAF")

#step 2: creation of sub-directories for the different GISTIC phases outputs
dir.create("1.rawdata_plots")
dir.create("2.Sep_plots")
dir.create("3.segment_data")
dir.create("4.segmented_plots")
dir.create("5.ASCAT_output_plots")


#_________PART 2: running ASCAT_________


#step 0: loading data
ascat.bc = ascat.loadData(file_logR, file_BAF, 
                          chrs = c(1:22))


#step 1: plotting raw logR and BAF data 
ascat.plotRawData(ascat.bc,
                  img.dir = "1.rawdata_plots")


#step 2: predict genotypes areas from BAF data
gg<-ascat.predictGermlineGenotypes(ascat.bc, 
                                   platform = "AffySNP6",
                                   img.dir = "2.Sep_plots")


#step 3: running ASPCF segmentation on both logR and BAF
ascat.bc = ascat.aspcf(ascat.bc, 
                       ascat.gg=gg, 
                       penalty=70, #default = 25 ( bigger => less segments)
                       out.dir= "3.segment_data") 


#step 4: plotting segmented data
ascat.plotSegmentedData(ascat.bc,
                        img.dir="4.segmented_plots")


#step 5: running ASCAT analysis and plotting outputs
ascat.output = ascat.runAscat(ascat.bc,
                              gamma=0.55, # default= 0.55
                              y_limit=5, # default= 5
                              circos = NA,
                              rho_manual = NA, # rho= aberrant cell fraction
                              psi_manual = NA, # psi= ploidy
                              img.dir = "5.ASCAT_output_plots")


#_________ part 3: exporting data results_________

# #visualize results
# ascat.output$aberrantcellfraction
# ascat.output$ploidy
# ascat.output$psi
# ascat.output$goodnessOfFit
# ascat.output$failedarrays
# ascat.output$nonaberrantarrays

# step 1: export a csv table of results
exportRes<- data.frame(ascat.output$goodnessOfFit,
                       ascat.output$aberrantcellfraction,
                       ascat.output$ploidy,
                       ascat.output$psi)

write.csv(exportRes, paste0( i, "_ASCAT.results_table.csv"))

#step 2: export the lists of failed arrays and non-aberrant arrays

write.csv(ascat.output$failedarrays, paste0( i,"_failed_arrays_list.csv"))
write.csv(ascat.output$nonaberrantarrays, paste0( i,"_nonAberrant_arrays_list.csv"))

#step 3: export a tsv file with the segments for every sample (adjusted-fitted and adjusted-raw)
library(readr)

segmentsAdj<-ascat.output$segments
segmentsAdj$CNvalue<- segmentsAdj$nMajor+segmentsAdj$nMinor
segmentsAdj$length<- segmentsAdj$endpos - segmentsAdj$startpos
write_tsv(segmentsAdj, paste0( i,"_adj_fitted_segments.tsv"))

segmentsAdjRaw<-ascat.output$segments_raw
segmentsAdjRaw$CNvalue<- segmentsAdjRaw$nMajor + segmentsAdjRaw$nMinor
segmentsAdjRaw$CNvalueRaw<- segmentsAdjRaw$nAraw + segmentsAdjRaw$nBraw
segmentsAdjRaw$length<- segmentsAdjRaw$endpos - segmentsAdjRaw$startpos
write_tsv(segmentsAdjRaw, paste0( i,"adj_fitted_and_raw_segments.tsv"))

# export ascat.output as a .Rdata single object - use readRDS() function to restore/load it in a new R session
saveRDS(ascat.output, paste0( i,"ascat_output.RData"))

# my_output<-readRDS("ascat_output.RData")





########################## BATCH 4: HD Failed #######################

# this batch is created from HD samples that failed the ASCAT run in batch 1 and 2: rerun them on SNP6 settings

library(ASCAT)

directories <- c("batchHDF")

i = directories[1]


#_________PART 1: set up dir tree and files_________

#step 1: choose working directory
# setwd(choose.dir()) 
setwd(paste0("D:/analisi_in_corso/risk_1q_13/ASCAT/", i))

#step 2: identify multi-sample logR and BAF files (pre-generated and already present in the working directory) 
file_logR<-list.files(pattern="logR")
file_BAF <-list.files(pattern="BAF")

#step 3: creation of sub-directories for the different GISTIC phases outputs
dir.create("1.rawdata_plots")
dir.create("2.Sep_plots")
dir.create("3.segment_data")
dir.create("4.segmented_plots")
dir.create("5.ASCAT_output_plots")


#_________PART 2: running ASCAT_________


#step 1: loading data
ascat.bc = ascat.loadData(file_logR, file_BAF)


#step 2: plotting raw logR and BAF data 
ascat.plotRawData(ascat.bc,
                  img.dir = "1.rawdata_plots")


#step 3: predict genotypes areas from BAF data
gg<-ascat.predictGermlineGenotypes(ascat.bc, 
                                   platform = "AffySNP6",
                                   img.dir = "2.Sep_plots")


#step 4: running ASPCF segmentation on both logR and BAF
ascat.bc = ascat.aspcf(ascat.bc, 
                       ascat.gg=gg, 
                       penalty=70, #default = 25 ( bigger => less segments)
                       out.dir= "3.segment_data") 


#step 5: plotting segmented data
ascat.plotSegmentedData(ascat.bc,
                        img.dir="4.segmented_plots")


#step 6: running ASCAT analysis and plotting outputs
ascat.output = ascat.runAscat(ascat.bc,
                              gamma=0.55, # default= 0.55
                              y_limit=5, # default= 5
                              circos = NA,
                              rho_manual = NA, # rho= aberrant cell fraction
                              psi_manual = NA, # psi= ploidy
                              img.dir = "5.ASCAT_output_plots")


#_________ part 3: exporting data results_________

# #visualize results
# ascat.output$aberrantcellfraction
# ascat.output$ploidy
# ascat.output$psi
# ascat.output$goodnessOfFit
# ascat.output$failedarrays
# ascat.output$nonaberrantarrays

#export a csv table of results
exportRes<- data.frame(ascat.output$goodnessOfFit,
                       ascat.output$aberrantcellfraction,
                       ascat.output$ploidy,
                       ascat.output$psi)

write.csv(exportRes, paste0( i, "_ASCAT.results_table.csv"))

#export the lists of failed arrays and non-aberrant arrays, if any

write.csv(ascat.output$failedarrays, paste0( i,"_failed_arrays_list.csv"))
write.csv(ascat.output$nonaberrantarrays, paste0( i,"_nonAberrant_arrays_list.csv"))

#export a .tsv with the segments for every sample (adjusted-fitted and adjusted-raw)
library(readr)

segmentsAdj<-ascat.output$segments
segmentsAdj$CNvalue<- segmentsAdj$nMajor+segmentsAdj$nMinor
segmentsAdj$length<- segmentsAdj$endpos - segmentsAdj$startpos
write_tsv(segmentsAdj, paste0( i,"_adj_fitted_segments.tsv"))

segmentsAdjRaw<-ascat.output$segments_raw
segmentsAdjRaw$CNvalue<- segmentsAdjRaw$nMajor + segmentsAdjRaw$nMinor
segmentsAdjRaw$CNvalueRaw<- segmentsAdjRaw$nAraw + segmentsAdjRaw$nBraw
segmentsAdjRaw$length<- segmentsAdjRaw$endpos - segmentsAdjRaw$startpos
write_tsv(segmentsAdjRaw, paste0( i,"adj_fitted_and_raw_segments.tsv"))

# export ascat.output as a .Rdata single object - use readRDS() function to restore/load it in a new R session
saveRDS(ascat.output, paste0( i,"ascat_output.RData"))

# my_outplut<-readRDS("ascat_output.RData")



########################## Experiment with GENDER ##########################################
# BATCH 1HD - 50 samples

library(ASCAT)
library(tidyverse)

#____ create the 50 samples files input from BATCH 1 _____

importLog <- data.table::fread("D:/analisi_in_corso/risk_1q_13/ASCAT/batchHD1/ascatReady_logR_1q13_batchHD1.txt")
importBAF <- data.table::fread("D:/analisi_in_corso/risk_1q_13/ASCAT/batchHD1/ascatReady_BAF_1q13_batchHD1.txt")

# cut the first 50 cols (+3 initial cols of annotations)
gender_50_log <- importLog[,1:53]
gender_50_BAF <- importBAF[,1:53]

# save the 50 sample files in a specific folder for a new ASCAT run 
readr::write_tsv(gender_50_log, "D:/analisi_in_corso/risk_1q_13/ASCAT/gender_run/ascatReady_logR_1q13_gender50.txt")
readr::write_tsv(gender_50_BAF, "D:/analisi_in_corso/risk_1q_13/ASCAT/gender_run/ascatReady_BAF_1q13_gender50.txt")

sampleNames <- names(gender_50_BAF)[-c(1:3)] %>% str_replace_all("^X", "") %>% str_replace_all("_\\.C", "_(C") %>% str_replace_all("y\\.", "y)") %>% str_replace_all("\\.", " ")
sampleNames

sexInfoFile <- readxl::read_xlsx("D:/analisi_in_corso/risk_1q_13/1q_13_clinicaldata_240519.xlsx") 

sexInfo <- sexInfoFile %>% select(nome_CEL_sample, Sex)

sexInfo50 <- sexInfo %>% filter(nome_CEL_sample %in% sampleNames)
sexInfo50 <- cbind(sexInfo50, sampleNames, names(gender_50_BAF)[-c(1:3)])

# # Checks
# sampleNames %in% sexInfo$nome_CEL_sample %>% table
# sampleNames[!(sampleNames %in% sexInfo$nome_CEL_sample)]


GENDER_VECTOR <- sexInfo50$Sex %>% str_replace("M", "XY") %>% str_replace("F", "XX")


directories <- c("gender_run")

i = directories[1]


#_________PART 1: set up dir tree and files_________

#step 1: choose working directory
# setwd(choose.dir()) 
setwd(paste0("D:/analisi_in_corso/risk_1q_13/ASCAT/", i))

#step 2: identify multi-sample logR and BAF files (pre-generated and already present in the working directory) 
file_logR<-list.files(pattern="logR")
file_BAF <-list.files(pattern="BAF")

#step 3: creation of sub-directories for the different GISTIC phases outputs
dir.create("1.rawdata_plots")
dir.create("2.Sep_plots")
dir.create("3.segment_data")
dir.create("4.segmented_plots")
dir.create("5.ASCAT_output_plots")


#_________PART 2: running ASCAT_________

#step 1: loading data
ascat.bc = ascat.loadData(file_logR, file_BAF, 
                          chrs = c(1:22,"X","Y","M"),
                          gender = GENDER_VECTOR) # ADDED the gender paramenter! 


#step 2: plotting raw logR and BAF data 
ascat.plotRawData(ascat.bc,
                  img.dir = "1.rawdata_plots")


#step 3: predict genotypes areas from BAF data
gg<-ascat.predictGermlineGenotypes(ascat.bc, 
                                   platform = "AffyCytoScanHD",
                                   img.dir = "2.Sep_plots")


#step 4: running ASPCF segmentation on both logR and BAF
ascat.bc = ascat.aspcf(ascat.bc, 
                       ascat.gg=gg, 
                       penalty=70, #default = 25 ( bigger => less segments)
                       out.dir= "3.segment_data") 


#step 5: plotting segmented data
ascat.plotSegmentedData(ascat.bc,
                        img.dir="4.segmented_plots")


#step 6: running ASCAT analysis and plotting outputs
ascat.output = ascat.runAscat(ascat.bc,
                              gamma=0.55, # default= 0.55
                              y_limit=5, # default= 5
                              circos = NA,
                              rho_manual = NA, # rho= aberrant cell fraction
                              psi_manual = NA, # psi= ploidy
                              img.dir = "5.ASCAT_output_plots")


#_________ part 3: exporting data results_________

# #visualize results
# ascat.output$aberrantcellfraction
# ascat.output$ploidy
# ascat.output$psi
# ascat.output$goodnessOfFit
# ascat.output$failedarrays
# ascat.output$nonaberrantarrays

#export a csv table of results
exportRes<- data.frame(ascat.output$goodnessOfFit,
                       ascat.output$aberrantcellfraction,
                       ascat.output$ploidy,
                       ascat.output$psi)

write.csv(exportRes, paste0( i, "_ASCAT.results_table.csv"))

#export the lists of failed arrays and non-aberrant arrays, if any

write.csv(ascat.output$failedarrays, paste0( i,"_failed_arrays_list.csv"))
write.csv(ascat.output$nonaberrantarrays, paste0( i,"_nonAberrant_arrays_list.csv"))

#export a .tsv with the segments for every sample (adjusted-fitted and adjusted-raw)
library(readr)

segmentsAdj<-ascat.output$segments
segmentsAdj$CNvalue<- segmentsAdj$nMajor+segmentsAdj$nMinor
segmentsAdj$length<- segmentsAdj$endpos - segmentsAdj$startpos
write_tsv(segmentsAdj, paste0( i,"_adj_fitted_segments.tsv"))

segmentsAdjRaw<-ascat.output$segments_raw
segmentsAdjRaw$CNvalue<- segmentsAdjRaw$nMajor + segmentsAdjRaw$nMinor
segmentsAdjRaw$CNvalueRaw<- segmentsAdjRaw$nAraw + segmentsAdjRaw$nBraw
segmentsAdjRaw$length<- segmentsAdjRaw$endpos - segmentsAdjRaw$startpos
write_tsv(segmentsAdjRaw, paste0( i,"adj_fitted_and_raw_segments.tsv"))

# export ascat.output as a .Rdata single object - use readRDS() function to restore/load it in a new R session
saveRDS(ascat.output, paste0( i,"ascat_output.RData"))



