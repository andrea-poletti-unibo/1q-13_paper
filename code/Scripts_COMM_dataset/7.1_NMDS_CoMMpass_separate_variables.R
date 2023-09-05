library(tidyverse)
library(data.table)
library(pca3d)
library(vegan)
library(ggrepel)
library(rgl)
library(gridExtra)


dir.create("plots/NMDS/CoMMpass/separate_variables/", showWarnings = F, recursive = T)

import<-fread("results/complete_database_CoMMpass_1q_13.txt")

# select only relevant variables to associate
names(import)

data<-import %>% select(Study_Visit_iD , contains("maj_broad"), contains("maj_focal"), contains("maj_call"), SeqWGS_WHSC1_CALL, SeqWGS_CCND1_CALL, SeqWGS_CCND3_CALL, SeqWGS_MAF_CALL, SeqWGS_MAFB_CALL)


names(data) %>% str_sort(numeric = T)


#____ KEEP ONLY UNIQUE/RELVANT CALLS____ 

# unique calls
data <- data %>% select( -c(AMP_maj_broad_chr_1q, 
                            AMP_maj_broad_chr_8q, 
                            DEL_maj_broad_chr_17p, 
                            DEL_maj_broad_chr_13q, 
                            DEL_maj_broad_chr_14q, 
                            DEL_maj_broad_chr_1p,
                            DEL_maj_broad_chr_16q, 
                            contains("maj_focal")
))

names(data)

# exclude calls from HyperDiploid chromosomes and single translocations for MDS dataset creation 

# renaming calls column
names(data) <- names(data) %>% str_replace("maj_broad_chr_|maj_call_chr_","")
names(data) <- names(data) %>% str_replace("SeqWGS","t")
names(data) <- names(data) %>% str_replace("_CALL","")
names(data)



##################################### NMDS ##################################### 

# metaMDS function

#________________ distances between alterations ____________

MDSdata <- t(as.matrix(data[,-1])) # transpose the data matrix (the method considers rows as the observation to localize), exclude sample column
colnames(MDSdata) <- data$sample

MDS.dist<- dist(MDSdata, method = "manhattan")

fit <- vegan::metaMDS(comm = MDS.dist, engine = "monoMDS", k=2, try= 50, trymax = 100)

ordiplot(fit, type = "text")
ordiplot(fit, type = "points")

ggdata <- as.data.frame(fit$points)
ggdata$alteration <- rownames(ggdata)

library(ggrepel)
ggplot(ggdata, aes(MDS1, MDS2)) + geom_point(size= 3, alpha=0.5) + 
  geom_label_repel( data = ggdata[ggdata$MDS2 < -50 | ggdata$MDS1>50 | ggdata$MDS1< -90  ,], aes(MDS1,MDS2, label=alteration)) 


ggplot(ggdata, aes(MDS1, MDS2)) + geom_point(size= 3, alpha=0.5) + 
  geom_text_repel( aes(MDS1,MDS2, label=alteration)) 


ggsave("plots/NMDS/CoMMpass/separate_variables/2D_NMDS_separatedVariables.png", height = 10, width = 12, units = "in")


#________________distances between samples __________________

MDSdata2 <- as.matrix(data[,-1]) # dont trasnpose, exclude sample column
rownames(MDSdata2) <- data$sample


MDS.dist2<- dist(MDSdata2, method = "manhattan")

fit2 <- vegan::metaMDS(comm = MDS.dist2, engine = "monoMDS", k=2, try= 50, trymax = 100)

ordiplot(fit2, type = "text")
ordiplot(fit2, type = "points")

category_HD <- ifelse(data$AMP_3p + data$AMP_3q + data$AMP_5p + data$AMP_5q + data$AMP_7p + data$AMP_7q + data$AMP_9p + data$AMP_9q + data$AMP_11p + data$AMP_11q+ data$AMP_15q+ data$AMP_19p+ data$AMP_19q+ data$AMP_21q > 3 ,"blue","gray10")
points(fit2, col=category_HD)

category_t_IgH <- ifelse(data$t_WHSC1 == 1 | data$t_CCND1 == 1 | data$t_CCND3 == 1 | data$t_MAF == 1 | data$t_MAFB == 1  ,"purple","gray10")
points(fit2, col=category_t_IgH)

category_AMP_1q <- ifelse(data$AMP_1q == 1 ,"green","gray10")
points(fit2, col=category_AMP_1q)

category_DEL_13q <- ifelse(data$DEL_13q == 1 ,"red","gray10")
points(fit2, col=category_DEL_13q)

risk_1q_13<- 3 - data$AMP_1q - data$DEL_13q

category_MMrisk <- ifelse(risk_1q_13 == 2 ,"gray60", ifelse( risk_1q_13 == 1, "orangered", "turquoise3"))
points(fit2, col=category_MMrisk)

category_risk1 <- ifelse(risk_1q_13 == 1 ,"orangered", "gray10")
points(fit2, col=category_risk1)

category_risk3 <- ifelse(risk_1q_13 == 3 ,"turquoise3", "gray10")
points(fit2, col=category_risk3)


# save these visualizations

#  NMDS ggplot 

FIT2 <- fit2$points %>% as.data.frame()

FIT2 %>% ggplot(aes(MDS1, MDS2)) + geom_point()

p1 <-FIT2 %>% ggplot(aes(MDS1, MDS2)) + geom_point( color=category_HD) + ggtitle("Hyperdiploidy - separated variables 2D NMDS")
p1

p2 <-FIT2 %>% ggplot(aes(MDS1, MDS2)) + geom_point( color=category_t_IgH) + ggtitle("t IgH - separated variables 2D NMDS")
p2

p3 <-FIT2 %>% ggplot(aes(MDS1, MDS2)) + geom_point( color=category_AMP_1q) + ggtitle("Amp 1q - separated variables 2D NMDS")
p3

p4 <-FIT2 %>% ggplot(aes(MDS1, MDS2)) + geom_point( color=category_DEL_13q) + ggtitle("Del 13 - separated variables 2D NMDS")
p4

p5 <-FIT2 %>% ggplot(aes(MDS1, MDS2)) + geom_point( color=category_MMrisk) + ggtitle("1q&13 classes - separated variables 2D NMDS")
p5

p6 <-FIT2 %>% ggplot(aes(MDS1, MDS2)) + geom_point( color=category_risk1) + ggtitle("1q&13+ vs others - separated variables 2D NMDS")
p6

# grid arrange

grid <- grid.arrange(p1,p2,p3,p4,p5,p6, ncol=3, nrow=2)
grid

ggsave(plot = grid, "plots/NMDS/CoMMpass/separate_variables/2D_NMDS_sepVar_grid_plot.png", height = 9, width = 14, units = "in")



#=============== 3d NMDS distances on samples ===================


fit3 <- vegan::metaMDS(comm = MDS.dist2, engine = "monoMDS", k=3, try= 50, trymax = 100)


# save.image("workfiles/NMDS_env_sepVar_CoMM.Rdata")
# load("workfiles/NMDS_env_sepVar_CoMM.Rdata") # for reproducibility


pca3d(fit3$points, col= category_HD, group= as.factor(category_HD), show.centroids = T , axe.titles=c("NMDS1","NMDS2","NMDS3"), shape="s" )

pca3d(fit3$points, col= category_t_IgH, group= as.factor(category_t_IgH), show.centroids = T , axe.titles=c("NMDS1","NMDS2","NMDS3"), shape="s" )

pca3d(fit3$points, col= category_MMrisk, group= as.factor(category_MMrisk), show.centroids = T , axe.titles=c("NMDS1","NMDS2","NMDS3"),  shape="s" )

pca3d(fit3$points, col= category_AMP_1q, group= as.factor(category_AMP_1q), show.centroids = T , axe.titles=c("NMDS1","NMDS2","NMDS3"),  shape="s" )

pca3d(fit3$points, col= category_DEL_13q, group= as.factor(category_DEL_13q), show.centroids = T , axe.titles=c("NMDS1","NMDS2","NMDS3"),  shape="s" )





########################## paper images #############################


outdir <- "plots/NMDS/CoMMpass/separate_variables/"
dir.create(outdir, showWarnings = F, recursive = T)


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
  
  snapshot3d(paste0(outdir,name,"_NMDS_CoMMpass_def3D.png"), height = 2000, width = 2000)

}




#____ center of Mass analysis _____

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



##################################################

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
          cex=1,
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

snapshot3d(paste0(outdir,"NMDS_CoMMpass_CLUSTERS.png"), height = 2000, width = 2000)
