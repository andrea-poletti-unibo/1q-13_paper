
############# from Commpass rawdata to GISTIC format #################

#__________ load libraries __________

library(data.table)
library(tidyverse)

#__________ data import e prep__________

segm <- fread("") # load Compass Segment data 

# filter per first visit (1_BM in sample name)

...
...
...

#__________ creation of a dataframe in GISTIC format___________

# Change the name of columns

GISTIC.segm<- data.frame( 
                       "Sample" = segm$___ , 
                       "Chromosome" = segm$___ , 
                       "Start Position" = segm$___ , 
                       "End Position" =segm$____ , 
                       "Num markers" = segm$_____ , 
                       "Seg.CN"= segm$____ 
                       )

# filter to exclude chr X, Y and Mitocondrial
GISTIC.segm.no_XYM<- filter(GISTIC.segm, Chromosome != "X", Chromosome != "Y", Chromosome != "M")


# save output
write_tsv(GISTIC.segm.no_XYM, "GISTIC_1q13_segm_CoMMpass.txt") 








