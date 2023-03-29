
library(tidyverse)
library(data.table)

df <- fread("results/complete_database_CoMMpass_1q_13.txt")

dfcalls <- df %>% select(matches("^DEL_|^AMP_"))


sums <- colSums(dfcalls)

res <- data.frame(event=names(dfcalls), n=sums, percentage=round(sums/nrow(df), 3))

res_maj <- res %>% filter(grepl("maj",res$event))
res_maj$event %>% as.character()

res_all <- res %>% filter(grepl("all",res$event))
res_all$event %>% as.character()

res_def <- cbind(res_all, res_maj)

names(res_def) <- c("event_all","n_all","percentage_all", "event_maj","n_maj", "percentage_maj")

event_minor <- res_def$event_all %>% str_replace("all","minor")
n_minor <- res_def$n_all - res_def$n_maj 
percentage_minor <- res_def$percentage_all - res_def$percentage_maj

res_def <- cbind(res_def, event_minor, n_minor, percentage_minor)

write_tsv(res_def,"results/CNA_calls_frequency_table_CoMM.txt")

