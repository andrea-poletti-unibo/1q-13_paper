library(data.table)
library(tidyverse)

db <- readxl::read_xlsx("D:/analisi_in_corso/risk_1q_13/1q_13_clinicaldata_240519.xlsx")

CELfiles <- db$nome_CEL_sample

Rawcopy <- list.dirs("D:/RAWCOPY/alpha= 10^-7/", recursive = F, full.names = F)

# CELfiles %in% Rawcopy %>% table


# create a separate list for snp HD and snp6.0 (different probeset length)

idx<- grep("CytoScanHD", CELfiles) # index for HD
filesHD<-CELfiles[idx]

idx<- grep("SNP_6", CELfiles) # index for 6.0
filesSNP6<-CELfiles[idx]


#========== run the extraction in BATCHES ==============

# LOOP to extract logR and BAF from the .Rdata rawcopy output of each sample
result.logR<-list()
result.BAF<-list()

# set base working directory
RawcopyDir <- "D:/RAWCOPY/alpha= 10^-7/"
basewd<-RawcopyDir

# DEFINE BATCHES (use 2 batches for HD samples and 1 for SNP6)

files1<-filesHD[1:190] # BATCH 1
files2<-filesHD[190:383] # BATCH 2
files3<-filesSNP6 # BATCH 3


for (i in files1) { # !!!define the batch!!!
  
  t1<-Sys.time()
  
  print(i)
  setwd(paste0(basewd,"/",i)) # change wd to a new sample
  load("rawcopy.Rdata") #load the sample specific .Rdata file
  result.logR[[i]]<-probes.txt$Value #add the sample logR to the total list
  result.BAF[[i]]<-probes.txt$B.Allele.Frequency # add the sample BAF to the total list
  
  t2<-Sys.time()
  print(t2-t1)
  
}

# transforming the list to a data frame
final.logR<-as.data.frame(result.logR) 
final.BAF<-as.data.frame(result.BAF)

# add columns with probe info
complete.logR<-base::cbind(probes.txt[,1:3],final.logR)
rm(final.logR)

complete.BAF<- base::cbind(probes.txt[,1:3],final.BAF)
rm(final.BAF)

# removing "chr" in Chromosome column
complete.logR$Chromosome<- gsub("chr","", complete.logR$Chromosome)
complete.BAF$Chromosome<- gsub("chr","", complete.BAF$Chromosome)

# renaming columns in ASCAT format
colnames(complete.logR)[1:3]<-c("probe","chrs","pos")
colnames(complete.BAF)[1:3]<-c("probe","chrs","pos")

# remove sex chromosomes and Mitochondrial 
LOGR_noXYM <- filter(complete.logR, chrs != "M",chrs != "X", chrs !="Y")
BAF_noXYM <- filter(complete.BAF, chrs != "M",chrs != "X", chrs !="Y")


# exporting files
setwd("D:/analisi_in_corso/risk_1q_13/ASCAT/batchHD1/")

library(readr)
write_tsv( LOGR_noXYM, "ascatReady.noXYM_logR_1q13_batchHD1.txt")
write_tsv( BAF_noXYM, "ascatReady.noXYM_BAF_1q13_batchHD1.txt")



####################### part 2: batch with failed HD samples ###############################

# this batch is created from HD samples that failed the ASCAT run in batch 1 and 2: rerun them on SNP6 settings

library(data.table)
library(tidyverse)

# import the files with the failed samples names from batch HD 1 and 2
HDF1 <- fread("D:/analisi_in_corso/risk_1q_13/ASCAT/batchHD1/batchHD1_failed_arrays_list.csv")[-1,] #remove first line (header)
HDF2 <- fread("D:/analisi_in_corso/risk_1q_13/ASCAT/batchHD2/batchHD2_failed_arrays_list.csv")[-1,] #remove first line (header)

HDF <- rbind(HDF1, HDF2)
# HDF$V2

# Reformat names to match rawcopy names
filesHDF <- HDF$V2 %>% str_replace("^X","") %>% str_replace("Array\\.","Array)") %>% str_replace("\\.","(")
# filesHDF

# check
# dirs <- list.dirs("D:/RAWCOPY/alpha= 10^-7/", recursive = F, full.names = F)
# filesHDF %in% dirs %>% table

# LOOP to extract logR and BAF from the .Rdata rawcopy output of each sample
result.logR<-list()
result.BAF<-list()

RawcopyDir <- "D:/RAWCOPY/alpha= 10^-7/"
basewd<-RawcopyDir

file<-filesHDF

for (i in file) {
  t1<-Sys.time()
  print(i)
  
  setwd(paste0(basewd,"/",i)) # change wd to a new sample
  load("rawcopy.Rdata") #load the sammple specific .Rdata file
  result.logR[[i]]<-probes.txt$Value #add the sample logR to the total list
  result.BAF[[i]]<-probes.txt$B.Allele.Frequency # add the sample BAF to the total list
  
  t2<-Sys.time()
  print(t2-t1)
}

final.logR<-as.data.frame(result.logR) # transforming the list to a data frame
final.BAF<-as.data.frame(result.BAF)

# add columns with probe info
complete.logR<-base::cbind(probes.txt[,1:3],final.logR)
rm(final.logR)

complete.BAF<- base::cbind(probes.txt[,1:3],final.BAF)
rm(final.BAF)

# removing "chr" in Chromosome column
complete.logR$Chromosome<- gsub("chr","", complete.logR$Chromosome)
complete.BAF$Chromosome<- gsub("chr","", complete.BAF$Chromosome)

# renaming columns in ASCAT format
colnames(complete.logR)[1:3]<-c("","chrs","pos")
colnames(complete.BAF)[1:3]<-c("","chrs","pos")

# exporting files
setwd("D:/analisi_in_corso/risk_1q_13/ASCAT/batchHDF/")

library(readr)
write_tsv(complete.logR,"ascatReady_logR_1q13_batchHDF.txt")
write_tsv(complete.BAF,"ascatReady_BAF_1q13_batchHDF.txt")
