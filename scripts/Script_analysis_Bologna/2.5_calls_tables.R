
library(tidyverse)

df <- data.table::fread("C:/Users/andre/Alma Mater Studiorum Università di Bologna/PROJECT 1q & 13 - Documenti/complete_database_1q_13_181019.txt")

dfcalls <- df %>% select(matches("^DEL|^AMP"))


sums <- colSums(dfcalls)

res <- data.frame(event=names(dfcalls), n=sums, percentage=round(sums/513, 3))

res_maj <- res %>% filter(grepl("maj",res$event))

res_all <- res %>% filter(grepl("all",res$event))

res_def <- cbind(res_all, res_maj)

names(res_def) <- c("event_all","n_all","percentage_all", "event_maj","n_maj", "percentage_maj")

event_minor <- res_def$event_all %>% str_replace("all","minor")
n_minor <- res_def$n_all - res_def$n_maj 
percentage_minor <- res_def$percentage_all - res_def$percentage_maj

res_def <- cbind(res_def, event_minor, n_minor, percentage_minor)

write.table(res_def, "clipboard", sep="\t", row.names=FALSE)

