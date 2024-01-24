
clean <- read.delim("ab24_count_TE.txt")
reference_table <- read.delim("reference_table.txt")
common_values <- merge(reference_table, clean, by = 'Name', all = TRUE) ##add NAs in mismatched
common_values1 <- merge(clean, reference_table, by = 'Name', all = TRUE)
common_values$Morc3ab_mut_Ave<-(c(common_values$morc3b.24hpf_rep1.T+common_values$morc3b.24hpf_rep2.T+common_values$morc3b.24hpf_rep3.T)/3)
common_values$Morc3ab_wt_Ave<-(c(common_values$wt_24hpf_rep1.C+common_values$wt_24hpf_rep2.C+common_values$wt_24hpf_rep3.C)/3)
write.table(common_values, file = "morc3b_24hpf_final.txt", sep = "\t", quote = F, row.names = T)
clean_final<-na.omit(common_values)
write.table(clean_final, "morcb24_final_TE.txt", sep = "\t", quote=F)
