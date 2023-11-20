
#####################  FISHER matrix SCRIPT ###########################

library(data.table)
library(tidyverse)

import<-fread("results/complete_database_CoMMpass_1q_13.txt")

rownames(import) <- import$Study_Visit_iD


# select only relevant variables to associate

data<-import %>% select(contains("maj_broad"), contains("maj_focal"), contains("maj_call"), HyperDiploidy, matches("SeqWGS"))
names(data)

#____ Build CALLS for traslocations ______

data$t_others <- data$SeqWGS_CCND2_CALL + data$SeqWGS_MAFB_CALL + data$SeqWGS_MAFA_CALL

#____ KEEP ONLY UNIQUE/RELVANT CALLS____ 

# unique calls
data <- data %>% select( -c(AMP_maj_broad_chr_1q, 
                            AMP_maj_broad_chr_8q, 
                            DEL_maj_broad_chr_17p, 
                            DEL_maj_broad_chr_13q, 
                            DEL_maj_broad_chr_14q, 
                            DEL_maj_broad_chr_1p,
                            DEL_maj_broad_chr_16q, 
                            contains("maj_focal"),
                            SeqWGS_MAFA_CALL,
                            SeqWGS_CCND2_CALL,
                            SeqWGS_MAFB_CALL
))
names(data)


data <- data %>% rename(t_WHSC1= SeqWGS_WHSC1_CALL,
                        t_CCND1= SeqWGS_CCND1_CALL,
                        t_MAF= SeqWGS_MAF_CALL,
                        t_MYC = SeqWGS_MYC_CALL)


names(data)

840*0.05

# relevant calls
data <- as.data.frame(data)
data %>% summarise_all(sum, na.rm =T) %>% as.numeric %>% sort %>% plot + abline(h = 10) + abline(h = 42, col="red")
rel <- data %>% summarise_all(sum, na.rm =T) %>% `>`(42) %>% as.vector

data <- data[rel]
names(data)
data %>% summarise_all(sum, na.rm =T) %>% `>`(42) %>% as.vector %>% table

 
#------------------ OPTIONAL: add group risk 1q & 13 ------------------------
MM_risk <- ifelse( data$AMP_maj_call_chr_1q == 1 & data$DEL_maj_call_chr_13q ==1, 1,
                   ifelse(data$AMP_maj_call_chr_1q == 1 | data$DEL_maj_call_chr_13q ==1, 2, 3) )

table(MM_risk)

data$MM_risk_1 <- ifelse(MM_risk==1, 1,0)
data$MM_risk_2 <- ifelse(MM_risk==2, 1,0)
data$MM_risk_3 <- ifelse(MM_risk==3, 1,0)

# dataprep COLUMNS

names(data) <- names(data) %>% str_replace("maj_broad_chr_|maj_call_chr_","") 


asd <- data

asd <- asd %>% select(str_order(colnames(asd), numeric = T))

# Fisher MATRIX
dataM <- as.matrix(asd)
rownames(dataM) <- rownames(import)
rownames(dataM)

sumsEvents <- apply(dataM, 2, sum, na.rm = TRUE)
sumsEvents

0.05*840


idxEvents <- which(sumsEvents >= 42)

dataMF <- dataM[,idxEvents]

apply(dataMF, 2, sum, na.rm = TRUE)

# table(dataMF[,1], dataMF[,3])
# f <- fisher.test(table(dataMF[,1], dataMF[,2]))
# c <- chisq.test(table(dataMF[,1], dataMF[,2]))
# 
# f$p.value
# c$p.value


i=1
j=3

mat <- matrix(data = NA, nrow= ncol(dataMF), ncol = ncol(dataMF), dimnames = list(colnames(dataMF),colnames(dataMF)))

for (i in 1:ncol(dataMF)){
  for(j in 1:ncol(dataMF)){
    f <- fisher.test(dataMF[,i], dataMF[,j], alternative = "greater")
    p <- f$p.value
    mat[i,j] <- p
  }
}


library(reshape2)

mat.UpTri<- mat
mat.UpTri[upper.tri(mat.UpTri, diag = T)] <- NA

mat.melt <- as.data.frame(melt(mat.UpTri))

mat.melt$logP <- log(mat.melt$value)

# create a thresholded version of p-value to handle the extremely significative values
mat.melt$pval.thresh <- mat.melt$value
mat.melt$pval.thresh[mat.melt$pval.thresh<0.001] <- 0.001

mat.melt$Var_fact <- mat.melt$Var1 %>% factor(levels=(mat.melt$Var1 %>% rev %>% unique))

ggplot(mat.melt, aes(mat.melt$Var_fact, mat.melt$Var2, fill=mat.melt$pval.thresh))+
  geom_tile(color="white")+
  scale_fill_gradient2(low="red", high="black", mid="orange", limits= c(0,0.05), midpoint = 0.025, space="Lab", name="Fisher p-value", na.value = "white")+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90, vjust=0,hjust=1))+
  coord_fixed()+
  labs(x="", y="") 



#=============  fisher exclusive assocition (inverted test) =================

mat_inv <- matrix(data = NA, nrow= ncol(dataMF), ncol = ncol(dataMF), dimnames = list(colnames(dataMF),colnames(dataMF)))
for (i in 1:ncol(dataMF)){
  for(j in 1:ncol(dataMF)){
    f <- fisher.test(dataMF[,i], dataMF[,j], alternative = "less")
    p <- f$p.value
    mat_inv[i,j] <- p
  }
}

mat.UpTri_inv<- mat_inv
mat.UpTri_inv[upper.tri(mat.UpTri_inv, diag = T)] <- NA

mat.melt_inv <- as.data.frame(melt(mat.UpTri_inv))

mat.melt_inv$logP <- log(mat.melt_inv$value)

mat.melt_inv$pval.thresh <- mat.melt_inv$value
mat.melt_inv$pval.thresh[mat.melt_inv$pval.thresh<0.001] <- 0.001


mat.melt_inv$Var_fact <- mat.melt_inv$Var1 %>% factor(levels=(mat.melt_inv$Var1 %>% rev %>% unique))

ggplot(mat.melt_inv, aes(mat.melt_inv$Var_fact, mat.melt_inv$Var2, fill=mat.melt_inv$pval.thresh))+
  geom_tile(color="white")+
  scale_fill_gradient2(low="blue", high="purple", mid="skyblue2", limits= c(0,0.05), midpoint = 0.025, space="Lab", name="Fisher p-value", na.value = "white")+
  theme_bw()+
  theme(axis.text.x=element_text(angle=60, vjust=1, size=10, hjust=1))+
  coord_fixed()+
  labs(x="", y="")




###################### ADD FUNCTIONS to add another scale ######################

#' Allows to add another scale 
#' 
#' @param new_aes character with the aesthetic for which new scales will be 
#' created



new_scale <- function(new_aes) {
  structure(ggplot2::standardise_aes_names(new_aes), class = "new_aes")
}

#' Convenient functions
new_scale_fill <- function() {
  new_scale("fill")
}

new_scale_color <- function() {
  new_scale("colour")
}

new_scale_colour <- function() {
  new_scale("colour")
}

#' Special behaviour of the "+" for adding a `new_aes` object
#' It changes the name of the aesthethic for the previous layers, appending
#' "_new" to them. 
ggplot_add.new_aes <- function(object, plot, object_name) {
  plot$layers <- lapply(plot$layers, bump_aes, new_aes = object)
  plot$scales$scales <- lapply(plot$scales$scales, bump_aes, new_aes = object)
  plot$labels <- bump_aes(plot$labels, new_aes = object)
  plot
}

bump_aes <- function(layer, new_aes) {
  UseMethod("bump_aes")
}

bump_aes.Scale <- function(layer, new_aes) {
  old_aes <- layer$aesthetics[remove_new(layer$aesthetics) %in% new_aes]
  new_aes <- paste0(old_aes, "_new")
  
  layer$aesthetics[layer$aesthetics %in% old_aes] <- new_aes
  
  if (is.character(layer$guide)) {
    layer$guide <- match.fun(paste("guide_", layer$guide, sep = ""))()
  }
  layer$guide$available_aes[layer$guide$available_aes %in% old_aes] <- new_aes
  layer
}

bump_aes.Layer <- function(layer, new_aes) {
  original_aes <- new_aes
  
  old_aes <- names(layer$mapping)[remove_new(names(layer$mapping)) %in% new_aes]
  new_aes <- paste0(old_aes, "_new")
  
  old_geom <- layer$geom
  
  old_setup <- old_geom$handle_na
  new_setup <- function(self, data, params) {
    colnames(data)[colnames(data) %in% new_aes] <- original_aes
    old_setup(data, params)
  }
  
  new_geom <- ggplot2::ggproto(paste0("New", class(old_geom)[1]), old_geom,
                               handle_na = new_setup)
  
  new_geom$default_aes <- change_name(new_geom$default_aes, old_aes, new_aes)
  new_geom$non_missing_aes <- change_name(new_geom$non_missing_aes, old_aes, new_aes)
  new_geom$required_aes <- change_name(new_geom$required_aes, old_aes, new_aes)
  new_geom$optional_aes <- change_name(new_geom$optional_aes, old_aes, new_aes)
  
  layer$geom <- new_geom
  
  old_stat <- layer$stat
  
  old_setup2 <- old_stat$handle_na
  new_setup <- function(self, data, params) {
    colnames(data)[colnames(data) %in% new_aes] <- original_aes
    old_setup2(data, params)
  }
  
  new_stat <- ggplot2::ggproto(paste0("New", class(old_stat)[1]), old_stat,
                               handle_na = new_setup)
  
  new_stat$default_aes <- change_name(new_stat$default_aes, old_aes, new_aes)
  new_stat$non_missing_aes <- change_name(new_stat$non_missing_aes, old_aes, new_aes)
  new_stat$required_aes <- change_name(new_stat$required_aes, old_aes, new_aes)
  new_stat$optional_aes <- change_name(new_stat$optional_aes, old_aes, new_aes)
  
  layer$stat <- new_stat
  
  layer$mapping <- change_name(layer$mapping, old_aes, new_aes)
  layer
}

bump_aes.list <- function(layer, new_aes) {
  old_aes <-  names(layer)[remove_new(names(layer)) %in% new_aes]
  new_aes <- paste0(old_aes, "_new")
  
  names(layer)[names(layer) %in% old_aes] <- new_aes
  layer
}

change_name <- function(list, old, new) {
  UseMethod("change_name")
}

change_name.character <- function(list, old, new) {
  list[list %in% old] <- new
  list
}

change_name.default <- function(list, old, new) {
  nam <- names(list)
  nam[nam %in% old] <- new
  names(list) <- nam
  list
}

change_name.NULL <- function(list, old, new) {
  NULL
}

remove_new <- function(aes) {
  stringi::stri_replace_all(aes, "", regex = "(_new)*")
}




#===== JOINT PLOT ========

mat.melt$pval.thresh_inv <- mat.melt_inv$pval.thresh

ggplot(mat.melt, aes(Var_fact, Var2))+
  theme_bw()+
  
  geom_tile(color="white", aes(fill=pval.thresh), data= mat.melt[mat.melt$pval.thresh<0.05 | is.na(mat.melt$pval.thresh_inv),])+
  scale_fill_gradient2(low="red", high="black", mid="orange", limits= c(0,0.05), midpoint = 0.025, space="Lab", name="Fisher p-value", na.value = "white")+
  
  new_scale("fill") +
  
  geom_tile(color="white", aes(fill= pval.thresh_inv),  data = mat.melt[mat.melt$pval.thresh_inv <0.05 | is.na(mat.melt$pval.thresh_inv) ,])+
  scale_fill_gradient2(low="blue", high="purple", mid="skyblue2", limits= c(0,0.05), midpoint = 0.025, space="Lab", name="Fisher p-value", na.value = "white")+
  
  theme(axis.text.x=element_text(angle=60, vjust=1, size=10, hjust=1),
        panel.border = element_blank(),
        legend.position = c(0.9, 0.8))+
  coord_fixed()+
  labs(x="", y="")



ggsave("plots/fisher_plots/fisher_matrix_CoMM.pdf", 
       dpi = 300, width = 10, height = 10, units = "in")

write_tsv(mat.melt, "plots/fisher_plots/fisher_matrix_CoMM_data.txt")
