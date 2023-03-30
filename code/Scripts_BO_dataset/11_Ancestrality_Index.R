library(data.table)
library(tidyverse)

import<-fread("results/complete_database_BO_1q_13.txt")

# select only relevant variables to associate
names(import)
data<-import %>% select(contains("maj_broad"), contains("maj_focal"), contains("maj_call"), matches("^t_"))
names(data) <- names(data) %>% str_replace("maj-", "maj_")
names(data)


#____ Build CALLS for traslocations ______

# data$t_others <- data$t_14_16 + data$t_14_20 + data$t_6_14

#____ KEEP ONLY UNIQUE/RELVANT CALLS____ 

data$IgH_traslocation <- ifelse(data$t_4_14 == 1 | data$t_6_14 ==1 | 
                                  data$t_11_14 ==1 | data$t_14_16== 1 | 
                                  data$t_14_20==1, 1,0)

data$IgH_traslocation %>% table

# unique calls
data <- data %>% select( -c(AMP_maj_broad_chr_1q, 
                            AMP_maj_broad_chr_8q, 
                            DEL_maj_broad_chr_17p, 
                            DEL_maj_broad_chr_13q, 
                            DEL_maj_broad_chr_14q, 
                            DEL_maj_broad_chr_1p,
                            DEL_maj_broad_chr_16q, 
                            contains("maj_focal"),
                            contains("t_")
))

names(data)

rownames(data) <- import$ID

dataC<-data %>% filter(complete.cases(data))


# compute the sum of alterations per patient (per row)
dataCS <- mutate(dataC, SUM=rowSums(dataC))


# inspect the samples with only one event
dataCS %>% filter(SUM==1) %>% pheatmap::pheatmap()

dataCS <- as.data.frame(dataCS)
rownames(dataCS) <- rownames(dataC)


dataCS <- rownames_to_column(dataCS, var = "sample")



#================ Loop to compute ancestrality index per sample ===============

timesMM<-list()
framesMM<-list()

i <- 2

# loop over every alteration 
for (i in 2:(ncol(dataCS)-1)) {
  
  message(names(dataCS)[i])
  
  # filter pts with the alteration
  p<-filter(dataCS, dataCS[,i] == 1)
  
  # compute mean total genomic complexity (MTGC) of those pts = mean of the sum of alterations per patient 
  MTGC<- mean(p$SUM)
  
  # save the number of those pts (frequency)
  number<- nrow(p)
  
  info<-c(MTGC, number)
  
  timesMM[[(names(dataCS)[i])]]<- info
  
  f<-data.frame("sample"= p$sample, "alteration"= rep(names(dataCS)[i],nrow(p)), "TGC"= p$SUM)
  
  framesMM[[i]]<-f
  
}

resultMM<-t(as.data.frame(timesMM))
colnames(resultMM)<-c("MTGC", "events")
resultMM <- as.data.frame(resultMM)

# compute penetrance / frequency of alterations
resultMM$penetrance<- round(resultMM$events/nrow(dataCS), 3)



#============ compute the Ancestrality Index ================

# AI = penetrance / MTGC

resultMM$Ancestrality_Index<- round(resultMM$penetrance / resultMM$MTGC, 5)* 100

resultMM <- rownames_to_column(resultMM,var = "alteration") 

exp <- resultMM %>% arrange(desc(Ancestrality_Index))

write_tsv(exp, "workfiles/AncestralityIndex_table_BO.txt")



