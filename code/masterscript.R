 setwd("..")
 
 # Load tidyverse
 library(tidyverse) # to pipe (%>%) and map across each file
 
 # List files and source each (BO dataset scripts folder)
 list.files("code/Scripts_BO_dataset/", full.names = TRUE) %>% map(source)
 
 # List files and source each (CoMMPass dataset scripts folder)
 list.files("code/Scripts_COMM_dataset/", full.names = TRUE) %>% map(source)