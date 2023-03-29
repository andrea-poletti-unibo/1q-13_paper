library(data.table)
library(tidyverse)


import <- fread("results/complete_database_CoMMpass_1q_13.txt")



# select only relevant variables to associate
names(import)
data<-import %>% select(contains("maj_broad"), contains("maj_focal"), matches("SeqWGS"), contains("maj_call"))
names(data)



#____ Build CALLS for traslocations ______

# data$t_others <- data$t_14_16 + data$t_14_20 + data$t_6_14

#____ KEEP ONLY UNIQUE/RELVANT CALLS____ 

data$IgH_traslocation <- ifelse(data$SeqWGS_WHSC1_CALL == 1 | data$SeqWGS_CCND1_CALL ==1 | 
                                  data$SeqWGS_CCND2_CALL ==1 | data$SeqWGS_CCND3_CALL== 1 | 
                                  data$SeqWGS_MAF_CALL==1 | data$SeqWGS_MAFA_CALL==1 | 
                                  data$SeqWGS_MAFB_CALL==1 | data$SeqWGS_WHSC1_CALL==1, 1,0)

data$IgH_traslocation %>% table


# unique calls
dataF <- data %>% select( -c(AMP_maj_broad_chr_1q, 
                            AMP_maj_broad_chr_8q, 
                            DEL_maj_broad_chr_17p, 
                            DEL_maj_broad_chr_13q, 
                            DEL_maj_broad_chr_14q, 
                            DEL_maj_broad_chr_1p,
                            DEL_maj_broad_chr_16q, 
                            contains("Seq"),
                            contains("maj_focal")
                            ))



names(dataF)

rownames(dataF) <- import$Study_Visit_iD

dataC<-dataF %>% filter(complete.cases(data))

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

write_tsv(resultMM, "workfiles/AncestralityIndex_table_CoMMpass.txt")


