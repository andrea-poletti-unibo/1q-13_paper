library(data.table)
library(tidyverse)


import<-fread("C:/Users/andre/Alma Mater Studiorum Università di Bologna/PROJECT 1q & 13 - Documenti/complete_database_1q_13_181019.txt")


# remove spaces from SNP names 
import$nSNP_1q_13_project
rownames(import) <- import$nSNP_1q_13_project %>% str_replace(" ", "_")
rownames(import)

import %>% names

# select only relevant variables to associate
names(import)
data<-import %>% select(contains("maj_broad"), contains("maj_focal"), matches("^t_[0-9]+"))
names(data) <- names(data) %>% str_replace("maj-", "maj_")
names(data)


#____ Build CALLS for arms with broad+focal_____

#amps
data$AMP_maj_call_chr_1q <- ifelse( data$AMP_maj_broad_chr_1q == 1 | data$AMP_maj_focal_ANP32E== 1 | data$AMP_maj_focal_CKS1B ==1 | data$AMP_maj_focal_MCL1 ==1, 1, 0)
data$AMP_maj_call_chr_8q <- ifelse( data$AMP_maj_broad_chr_8q == 1 | data$AMP_maj_focal_MYC == 1, 1, 0)
#dels
data$DEL_maj_call_chr_17p <- ifelse( data$DEL_maj_broad_chr_17p == 1 | data$DEL_maj_focal_TP53 == 1, 1, 0)
data$DEL_maj_call_chr_1p <-  ifelse( data$DEL_maj_broad_chr_1p == 1 | data$DEL_maj_focal_CDKN2C == 1 | data$DEL_maj_focal_FAM46C == 1 | data$DEL_maj_focal_FAF1 == 1,  1, 0)
data$DEL_maj_call_chr_13q <- ifelse( data$DEL_maj_broad_chr_13q == 1 | data$DEL_maj_focal_RB1 == 1, 1, 0)
data$DEL_maj_call_chr_14q <- ifelse( data$DEL_maj_broad_chr_14q == 1 | data$DEL_maj_focal_TRAF3 == 1, 1, 0)
data$DEL_maj_call_chr_16q <- ifelse( data$DEL_maj_broad_chr_16q == 1 | data$DEL_maj_focal_CYLD == 1, 1, 0)

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

rownames(data) <- rownames(import)

dataC<-data %>% filter(complete.cases(data))

rownames(dataC)

dataCS <- mutate(dataC, SUM=rowSums(dataC))

dataCS %>% filter(SUM==1) %>% pheatmap::pheatmap()

dataCS <- as.data.frame(dataCS)

rownames(dataCS) <- rownames(dataC)

dataCS <- rownames_to_column(dataCS, var = "sample")

timesMM<-list()
framesMM<-list()

i <- 1
for (i in 2:(ncol(dataCS)-1)) {
  print(i)
  p<-filter(dataCS, dataCS[,i] == 1)
  TO<- mean(p$SUM)
  number<- nrow(p)
  info<-c(TO, number)
  timesMM[[(names(dataCS)[i])]]<- info
  
  f<-data.frame("sample"= p$sample, "alteration"= rep(names(dataCS)[i],nrow(p)), "TGC"= p$SUM)
  framesMM[[i]]<-f
  
}



# risultati delle mediane dei valori di complessità genomica totale per ogni alterazione
resultMM<-t(as.data.frame(timesMM))
colnames(resultMM)<-c("MTGC", "events")
resultMM<-as.data.frame(resultMM)

resultMM$penetrance<- round(resultMM$events/nrow(dataCS), 3)

resultMM$Ancestrality_Index<- round(resultMM$penetrance / resultMM$MTGC, 5)* 100

resultMM <- rownames_to_column(resultMM,var = "alteration")

write_tsv(resultMM, "AncestralityIndex_table_Bologna_060821.txt")


######################################################################################


#creazione data frame con i singoli valori di MTGC per plot
resultF.MM<-Reduce(rbind, framesMM)

resultF.MM$alteration %>% as.character() %>% unique() %>% length()

library(ggplot2)
library(ggridges)
#plot histogram - MTCG per pts
ggplot(resultF.MM, aes(x = TGC, group= as.factor(alteration))) +
  geom_histogram(fill="blue", binwidth=1) +
  # scale_fill_cyclical(values = c("green3", "palegreen3")) +
  facet_wrap(alteration~., ncol = 3, nrow = 26)+
  xlim(c(0,35))

ggsave(filename = "C:/Users/andre/Alma Mater Studiorum Università di Bologna/PROJECT 1q & 13 - Documenti/Supplementary/Ancestrality_index_bolo_all_events_060821.png", device = "png", width = 12, height = 36)




