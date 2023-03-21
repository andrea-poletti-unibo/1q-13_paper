library(tidyverse)
library(data.table)
library(rgl)
library(pca3d)
library(MASS)
library(vegan)
library(ggrepel)
library(plot3D)

# setwd("D:/analisi_in_corso/risk_1q_13/")
# import <- fread("complete_database_1q_13_181019.txt")

import<-fread("C:/Users/andre/Alma Mater Studiorum Università di Bologna/PROJECT 1q & 13 - Documenti/complete_database_1q_13_181019.txt")

rownames(import) <- import$nSNP_1q_13_project %>% str_replace(" ", "_")
import <- rownames_to_column(import, "sample")

# select only relevant variables to associate
names(import)
data<-import %>% dplyr::select(sample, contains("maj_broad"), contains("maj_focal"), HyperDiploidy, matches("^t_"))
names(data)

#____ Build CALLS for arms with broad+focal_____

#amps
data$AMP_maj_call_chr_1q <- ifelse( data$AMP_maj_broad_chr_1q == 1 | data$AMP_maj_focal_ANP32E== 1 | data$AMP_maj_focal_CKS1B ==1 | data$AMP_maj_focal_MCL1 ==1, 1, 0)
data$AMP_maj_call_chr_8q <- ifelse( data$AMP_maj_broad_chr_8q == 1 | data$AMP_maj_focal_MYC == 1, 1, 0)

#dels
data$DEL_maj_call_chr_17p <- ifelse( data$DEL_maj_broad_chr_17p == 1 | data$DEL_maj_focal_TP53 == 1, 1, 0)
data$DEL_maj_call_chr_1p <- ifelse( data$DEL_maj_broad_chr_1p == 1 | data$DEL_maj_focal_CDKN2C == 1 | data$DEL_maj_focal_FAM46C == 1 | data$DEL_maj_focal_FAF1 == 1,  1, 0)
data$DEL_maj_call_chr_13q <- ifelse( data$DEL_maj_broad_chr_13q == 1 | data$DEL_maj_focal_RB1 == 1, 1, 0)
data$DEL_maj_call_chr_14q <- ifelse( data$DEL_maj_broad_chr_14q == 1 | data$DEL_maj_focal_TRAF3 == 1, 1, 0)
data$DEL_maj_call_chr_16q <- ifelse( data$DEL_maj_broad_chr_16q == 1 | data$DEL_maj_focal_CYLD == 1, 1, 0)

#____ Build CALLS for traslocations ______

data$t_others <- data$t_14_16 + data$t_14_20 + data$t_6_14

#____ KEEP ONLY UNIQUE/RELVANT CALLS____ 

# unique calls
data <- data %>% dplyr::select( -c(AMP_maj_broad_chr_1q, 
                            AMP_maj_broad_chr_8q, 
                            DEL_maj_broad_chr_17p, 
                            DEL_maj_broad_chr_13q, 
                            DEL_maj_broad_chr_14q, 
                            DEL_maj_broad_chr_1p,
                            DEL_maj_broad_chr_16q, 
                            contains("maj_focal"),
                            t_14_16,
                            t_6_14,
                            t_14_20
))

names(data)

# exclude calls from HyperDiploid chromosomes and single translocations for MDS dataset creation 
data <- data %>% dplyr::select( -c(AMP_maj_broad_chr_3p, 
                            AMP_maj_broad_chr_3q, 
                            AMP_maj_broad_chr_5p, 
                            AMP_maj_broad_chr_5q, 
                            AMP_maj_broad_chr_7p, 
                            AMP_maj_broad_chr_7q,
                            AMP_maj_broad_chr_9p, 
                            AMP_maj_broad_chr_9q,
                            AMP_maj_broad_chr_11p, 
                            AMP_maj_broad_chr_11q,
                            AMP_maj_broad_chr_15q,
                            AMP_maj_broad_chr_19p,
                            AMP_maj_broad_chr_19q,
                            AMP_maj_broad_chr_21q,
                            t_4_14,
                            t_11_14,
                            t_others
                            ))

# renaming calls column
names(data) <- names(data) %>% str_replace("maj_broad_chr_|maj_call_chr_","")
names(data)














MDSdata <- data

# filter MDSdata for missing data (missing traslocation from AMB or Bo2005 protocols)

missing_T <- is.na(import$t_11_14) + 
  is.na(import$t_4_14) + 
  is.na(import$t_6_14) + 
  is.na(import$t_14_16) + 
  is.na(import$t_14_20)
table(missing_T)
import$missing_T <- missing_T

table(missing_T, import$t_IgH, exclude = NULL)
table(missing_T, import$PROTOCOL, exclude = NULL)
table(import$t_IgH, import$PROTOCOL, exclude = NULL)

# import %>% dplyr::select(missing_T, t_IgH, PROTOCOL) %>% View

df.2<- dplyr::filter(MDSdata, !is.na(t_IgH))


data$HyperDiploidy %>% table
295/513



library(MASS)
# swiss.x <- as.matrix(swiss[, -1])
# swiss.dist <- dist(swiss.x)
# 
# swiss.mds <- isoMDS(swiss.dist)
# head(swiss.mds$points)
# 
# plot(swiss.mds$points, type = "n")
# text(swiss.mds$points, labels = as.character(1:nrow(swiss.x)))
# 
# swiss.sh <- Shepard(swiss.dist, swiss.mds$points)
# plot(swiss.sh, pch = ".")
# 
# lines(swiss.sh$x, swiss.sh$yf, type = "S")
####################################

# isoMDS function

# MDSdata <- as.matrix(df.2)
# 
# MDS.dist<- dist(t(MDSdata), method = "manhattan")
# 
# str(MDS.dist)
# summary(MDS.dist)
# 
# MDS<- isoMDS(MDS.dist)

####################################

# metaMDS function
library(vegan)

#________________ distances between alterations ____________
MDSdata <- t(as.matrix(df.2[,-1])) # transpose the data matrix (the method considers rows as the observation to localize), exclude sample column
colnames(MDSdata) <- df.2$sample

MDS.dist<- dist(MDSdata, method = "manhattan")

fit <- vegan::metaMDS(comm = MDS.dist, engine = "monoMDS", k=2, try= 50, trymax = 100)

ordiplot(fit, type = "text")
ordiplot(fit, type = "points")

ggdata <- as.data.frame(fit$points)
ggdata$alteration <- rownames(ggdata)

library(ggrepel)
ggplot(ggdata, aes(MDS1, MDS2)) + geom_point(size= 3, alpha=0.5) + 
  geom_label_repel( data = ggdata[ggdata$MDS1 > 50 | ggdata$MDS1 < -50,], aes(MDS1,MDS2, label=alteration)) 

ggplot(ggdata, aes(MDS1, MDS2)) + geom_point(size= 3, alpha=0.5) + 
  geom_label_repel( aes(MDS1,MDS2, label=alteration)) 



#________________distances between samples __________________

MDSdata2 <- as.matrix(df.2[,-1]) # dont trasnpose, exclude sample column
rownames(MDSdata2) <- df.2$sample


MDS.dist2<- dist(MDSdata2, method = "manhattan")

fit2 <- vegan::metaMDS(comm = MDS.dist2, engine = "monoMDS", k=2, try= 50, trymax = 100)

ordiplot(fit2, type = "text")
ordiplot(fit2, type = "points")

category_HD <- ifelse(df.2$HyperDiploidy == 1 ,"blue","gray10")
points(fit2, col=category_HD)

category_t_IgH <- ifelse(df.2$t_IgH == 1 ,"purple","gray10")
points(fit2, col=category_t_IgH)

category_AMP_1q <- ifelse(df.2$AMP_1q == 1 ,"green","gray10")
points(fit2, col=category_AMP_1q)

category_DEL_13q <- ifelse(df.2$DEL_13q == 1 ,"red","gray10")
points(fit2, col=category_DEL_13q)

risk_1q_13<- 3 - df.2$AMP_1q - df.2$DEL_13q

category_MMrisk <- ifelse(risk_1q_13 == 2 ,"gray60", ifelse( risk_1q_13 == 1, "orangered", "turquoise3"))
points(fit2, col=category_MMrisk)

category_risk1 <- ifelse(risk_1q_13 == 1 ,"green", "red")
points(fit2, col=category_risk1)

category_risk3 <- ifelse(risk_1q_13 == 3 ,"turquoise3", "red")
points(fit2, col=category_risk3)

# NMDS ggplot 

FIT2 <- fit2$points %>% as.data.frame()

FIT2 %>% ggplot(aes(MDS1, MDS2)) + geom_point()




#=============== 3d NMDS distances on samples ===================
library(rgl)
library(pca3d)
library(vegan3d)

fit3 <- vegan::metaMDS(comm = MDS.dist2, engine = "monoMDS", k=3, try= 50, trymax = 100)

fit3$grstress # stress: 0.1251906

# save.image("Script/4.2_NMDS_env.Rdata")
load("Script/4.2_NMDS_env.Rdata")


# shape <-  ifelse(risk_1q_13 == 1 ,"t", ifelse( risk_1q_13 == 2, "s", "c"))

outdir <- "C:/Users/andre/Alma Mater Studiorum Università di Bologna/PROJECT 1q & 13 - Documenti/PAPER_FIGURES_OFFICIAL/"

#_________________ paper images ________________________

pca3d(fit3$points, col= category_HD, group= as.factor(category_HD), 
      show.centroids = T , axe.titles=c("NMDS1","NMDS2","NMDS3"), shape="s")
# snapshotPCA3d(paste0(outdir,"NMDS_bolo_hyperdiploid.png"))
# rgl.snapshot("C:/Users/andre/Desktop/test.png", width=3000, height=3000)

snapshot3d(paste0(outdir,"NMDS_bolo_hyperdiploid.png"), width=10000, height=10000)


rgl.postscript("C:/Users/mm_gr/Desktop/test.pdf",fmt="pdf")


view3d( theta = 0, phi = 15)

pca3d(fit3$points, col= category_t_IgH, group= as.factor(category_t_IgH), 
      show.centroids = T , axe.titles=c("NMDS1","NMDS2","NMDS3"), shape="s" )
snapshot3d(paste0(outdir,"NMDS_bolo_t_IGH.png"), width=10000, height=10000)


pca3d(fit3$points, col= category_MMrisk, group= as.factor(category_MMrisk), 
      show.centroids = T , axe.titles=c("NMDS1","NMDS2","NMDS3"),  shape="s" )
snapshot3d(paste0(outdir,"NMDS_bolo_MMgroups.png"), width=10000, height=10000)


pca3d(fit3$points, col= category_AMP_1q, group= as.factor(category_AMP_1q), 
      show.centroids = T , axe.titles=c("NMDS1","NMDS2","NMDS3"),  shape="s" )
snapshot3d(paste0(outdir,"NMDS_bolo_amp1q.png"), width=10000, height=10000)


pca3d(fit3$points, col= category_DEL_13q, group= as.factor(category_DEL_13q), 
      show.centroids = T , axe.titles=c("NMDS1","NMDS2","NMDS3"),  shape="s" )
snapshot3d(paste0(outdir,"NMDS_bolo_del13.png"), width=10000, height=10000)



#________________________________________________________
# other

risk_1q_13<- df.2$AMP_1q + df.2$DEL_13q
category <- ifelse(risk_1q_13 == 1 ,"gray80", ifelse( risk_1q_13 == 2, "orangered", "turquoise3"))
pca3d(fit3$points, col= category, group= as.factor(category), show.centroids = T, axe.titles=c("NMDS1","NMDS2","NMDS3"), shape="s" )

category <- ifelse(risk_1q_13 == 2 ,"green", "red")
pca3d(fit3$points, col= category, group= as.factor(category), show.centroids = T, axe.titles=c("NMDS1","NMDS2","NMDS3"), shape="s" )

category <- ifelse(risk_1q_13 == 1 ,"turquoise3", "red")
pca3d(fit3$points, col= category, group= as.factor(category), show.centroids = T, axe.titles=c("NMDS1","NMDS2","NMDS3"), shape="s" )

category <- ifelse(risk_1q_13 == 0 ,"turquoise3", "red")
pca3d(fit3$points, col= category, group= as.factor(category), show.centroids = T, axe.titles=c("NMDS1","NMDS2","NMDS3"), shape="s" )


############################ paper images PT2 #############################

library(rgl)
load("Script/4.2_NMDS_env.Rdata")


category_risk1 <- ifelse(risk_1q_13 == 1 ,"orangered", "gray10")

category_risk3 <- ifelse(risk_1q_13 == 3 ,"turquoise3", "gray10")

CATS <- list(
  HD=category_HD,
  t_IgH=category_t_IgH,
  AMP_1q= category_AMP_1q,
  DEL_13=category_DEL_13q,
  MM_risk=category_MMrisk,
  risk_1=category_risk1,
  risk_3=category_risk3
  )

COLS <- list(
  HD=c("blue","grey20"),
  t_IgH=c("grey20","purple"),
  AMP_1q= c("grey20","green"),
  DEL_13=c("grey20","red"),
  MM_risk=c("gray50","orangered","turquoise3"),
  risk_1=c("grey20","orangered"),
  risk_3=c("grey20","turquoise3")
  )


i=1

for(i in seq_along(CATS)){
  
  i
  name <- names(CATS[i])
  print(name)
  
  categ <- CATS[i] %>% unlist()
  groups <- CATS[i] %>% unlist %>%  as.factor()
  levs <- levels(groups)
  group.col <- COLS[name] %>% unlist
  ef <- 1.3
  
  df <- fit3$points
  x <- df[,1]
  y <- df[,2]
  z <- df[,3]
  
  rgl.open() # Open a new RGL device
  
  rgl.bg(color = "white") # Setup the background color
  
  rgl.spheres(df, r = 0.25, color = categ) # ADD POINTS
  
  lim <- function(x){c(min(x), max(x)) * ef}
  # Add axes
  rgl.lines(lim(x), c(0, 0), c(0, 0), color = "grey70", lwd=3)
  rgl.lines(c(0, 0), lim(y), c(0, 0), color = "grey70", lwd=3)
  rgl.lines(c(0, 0), c(0, 0), lim(z), color = "grey70", lwd=3)
  
  # Add a point at the end of each axes to specify the direction
  axes <- rbind(c(lim(x)[2], 0, 0), 
                c(0, lim(y)[2], 0), 
                c(0, 0, lim(z)[2]))
  rgl.points(axes, size=9, color = "grey70") 
  
  # Add axis labels
  rgl.texts(axes, text = c("NMDS1", "NMDS2", "NMDS3"), 
            color = "grey70",
            cex=2,
            adj = c(0.5, -0.8), 
            size = 2)
  
  # # Add plane
  # xlim <- lim(x)/ef
  # zlim <- lim(z)/ef
  # rgl.quads( x = rep(xlim, each = 2),
  #            y = c(0, 0, 0, 0),
  #            z = c(zlim[1], zlim[2], zlim[2], zlim[1]),
  #            alpha=0.15)
  
  # # Compute ellipse for each group
  # for (j in 1:length(levs)) {
  #   group <- levs[j]
  #   selected <- groups == group
  #   xx <- x[selected]
  #   yy <- y[selected]
  #   zz <- z[selected]
  #   ellips <- ellipse3d(cov(cbind(xx,yy,zz)),
  #                       centre=c(mean(xx), mean(yy), mean(zz)), level = 0.95)
  #   shade3d(ellips, col = group.col[j], alpha = 0.1, lit = F)
  #   wire3d(ellips, col = group.col[j], alpha = 0.3, lit = F)
  # }
  
  # Compute centers and lines for each group
  for (j in 1:length(levs)) {
    group <- levs[j]
    selected <- groups == group
    xx <- x[selected] 
    yy <- y[selected]
    zz <- z[selected]
    
    X1 <- mapply(c, xx, mean(xx), SIMPLIFY = F) %>% unlist
    Y1 <- mapply(c, yy, mean(yy), SIMPLIFY = F) %>% unlist
    Z1 <- mapply(c, zz, mean(zz), SIMPLIFY = F) %>% unlist
    
    rgl.lines(x=X1, y=Y1, z=Z1, color = group.col[j], size=1)
  }
  
  # good POV
  rgl.viewpoint(theta = 60, phi = 30, fov = 60, zoom = 0.7)
  
  # standard POV
  # rgl.viewpoint(theta = 45, phi = 45, fov = 60, zoom = 0.7)
  
  snapshot3d(paste0("C:/Users/andre/Alma Mater Studiorum Università di Bologna/PROJECT 1q & 13 - Documenti/PAPER_FIGURES_OFFICIAL/NMDS/",name,"_NMDS_Bolo_def3D.png"), height = 2000, width = 2000)
}



#____ center of Mass analysis _____



name <- names(CATS[i])
print(name)
categ <- CATS[i] %>% unlist()
groups <- CATS[i] %>% unlist %>%  as.factor()
levs <- levels(groups)
group.col <- COLS[name] %>% unlist


ef <- 1.3
df <- fit3$points
x <- df[,1]
y <- df[,2]
z <- df[,3]

rgl.open() # Open a new RGL device

rgl.bg(color = "white") # Setup the background color

lim <- function(x){c(min(x), max(x)) * ef}
# Add axes
rgl.lines(lim(x), c(0, 0), c(0, 0), color = "grey70", lwd=3)
rgl.lines(c(0, 0), lim(y), c(0, 0), color = "grey70", lwd=3)
rgl.lines(c(0, 0), c(0, 0), lim(z), color = "grey70", lwd=3)

# Add a point at the end of each axes to specify the direction
axes <- rbind(c(lim(x)[2], 0, 0), 
              c(0, lim(y)[2], 0), 
              c(0, 0, lim(z)[2]))
rgl.points(axes, size=9, color = "grey70") 

# Add axis labels
rgl.texts(axes, text = c("NMDS1", "NMDS2", "NMDS3"), 
          color = "grey70",
          cex=2,
          adj = c(0.5, -0.8), 
          size = 2)



for(i in seq_along(CATS)){
  
  name <- names(CATS[i])
  print(name)
  categ <- CATS[i] %>% unlist()
  groups <- CATS[i] %>% unlist %>%  as.factor()
  levs <- levels(groups)
  group.col <- COLS[name] %>% unlist
  
  # Compute centers and lines for each group
  for (j in 1:length(levs)) {
    group <- levs[j]
    selected <- groups == group
    xx <- x[selected] 
    yy <- y[selected]
    zz <- z[selected]
    
    rgl.spheres(mean(xx), mean(yy), mean(zz), color=levs[j], r=0.5)
    
  }
}
# good POV
rgl.viewpoint(theta = 60, phi = 30, fov = 60, zoom = 0.7)



######################################


i=1

ef <- 1.3
df <- fit3$points
x <- df[,1]
y <- df[,2]
z <- df[,3]


res <- data.frame()
for(i in seq_along(CATS)){
  name <- names(CATS[i])
  print(name)
  categ <- CATS[i] %>% unlist()
  groups <- CATS[i] %>% unlist %>%  as.factor()
  levs <- levels(groups)
  group.col <- COLS[name] %>% unlist
  
  # Compute centers and lines for each group
  j=1
  for (j in 1:length(levs)) {
    group <- levs[j]
    selected <- groups == group
    xx <- x[selected] 
    yy <- y[selected]
    zz <- z[selected]
    
    df <- data.frame(name=name, group=group, x=mean(xx), y= mean(yy), z=mean(zz))
    res <- rbind(res,df)
  }
}




res2 <- res %>% dplyr::filter(group!="gray10", name %in% c("HD", "t_IgH", "risk_1", "risk_3") )
rownames(res2) <- res2$name

rgl.open() # Open a new RGL device
rgl.bg(color = "white") # Setup the background color

x=res2$x
y=res2$y
z=res2$z
rgl.spheres(x, y, z, r=0.1, color = res2$group) 

lim <- function(x){c(min(x), max(x)) * ef}
# Add axes
rgl.lines(lim(x), c(0, 0), c(0, 0), color = "grey70", lwd=3)
rgl.lines(c(0, 0), lim(y), c(0, 0), color = "grey70", lwd=3)
rgl.lines(c(0, 0), c(0, 0), lim(z), color = "grey70", lwd=3)

# Add a point at the end of each axes to specify the direction
axes <- rbind(c(lim(x)[2], 0, 0), 
              c(0, lim(y)[2], 0), 
              c(0, 0, lim(z)[2]))
rgl.points(axes, size=9, color = "grey70") 

# Add axis labels
rgl.texts(axes, text = c("NMDS1", "NMDS2", "NMDS3"), 
          color = "grey70",
          cex=2,
          adj = c(0.5, -0.8), 
          size = 2)

res2["HD","x"]

rgl.lines(c(res2["HD","x"], res2["t_IgH","x"]), 
          c(res2["HD","y"],res2["t_IgH","y"]),
          c(res2["HD","z"], res2["t_IgH","z"]), 
          color = "grey20", lwd=9)

rgl.lines(c(res2["risk_1","x"], res2["risk_3","x"]), 
          c(res2["risk_1","y"],res2["risk_3","y"]),
          c(res2["risk_1","z"], res2["risk_3","z"]), 
          color = "grey20", lwd=9)

# rgl.lines(c(res2["HD","x"], res2["HD","x"]), 
#           c(res2["HD","y"], 0),
#           c(res2["HD","z"], res2["HD","z"]), 
#           color = "grey80", lwd=4)

# Add plane
xlim <- lim(x)/1.2
zlim <- lim(z)/1.2
rgl.quads( x = rep(xlim, each = 2),
           y = c(0, 0, 0, 0),
           z = c(zlim[1], zlim[2], zlim[2], zlim[1]),
           alpha=0.2)

rgl.viewpoint(theta = 60, phi = 30, fov = 60, zoom = 0.7)

snapshot3d(paste0("C:/Users/andre/Alma Mater Studiorum Università di Bologna/PROJECT 1q & 13 - Documenti/PAPER_FIGURES_OFFICIAL/NMDS/",name,"_NMDS_Bolo_CLUSTERS.png"), height = 2000, width = 2000)
