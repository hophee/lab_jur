suppressPackageStartupMessages(library(tidyverse))
library(stringi)

arg <- commandArgs(trailingOnly = TRUE)

meta <- read_tsv('metadata.tsv', show_col_types = FALSE, name_repair = "universal_quiet") %>% janitor::clean_names()
meta <- meta %>% select(1, 2, 5, 6, 7, 8) %>% group_by(familia)
meta$familia <- iconv(stri_trans_general(meta$familia, "russian-latin/bgn"),from="UTF8",to="ASCII",sub="") %>% tolower()

out_csv <- data.frame(patient=(rep(unique(meta$familia), each =2))) %>% 
  mutate(sex='XX', status=rep(0:1, 28), sample=paste0(rep(c('normal_', 'tumor_'), 28), patient), lane=paste0('lane_', 1:56))
orig_fams <- meta %>% pull(familia) %>% unique()
ids_norm <- meta %>% filter(nukleinovaa_kislota=='Герминальная ДНК') %>% select(familia, id_obrazca)%>% rename('norm_id'='id_obrazca')
ids_tumor <- meta %>% filter(nukleinovaa_kislota=='Соматическая ДНК') %>% select(familia, id_obrazca) %>% rename('tumor_id'='id_obrazca')
ids_tab <- left_join(ids_norm, ids_tumor, by='familia') %>% rename('patient'='familia')
ids <- vector(length = 56)
ids[seq(1, 56, 2)] <- ids_tab$norm_id
ids[seq(2, 56, 2)] <- ids_tab$tumor_id
out_csv <-  out_csv %>% mutate(fastq_1=paste0(c('exome_norm/', 'exome_tumor/'),ids, '_1.fq.gz'), fastq_2=paste0(c('exome_norm/', 'exome_tumor/'),ids, '_2.fq.gz'))

arg_num <- as.numeric(arg)
iwgr <- arg_num[1]:arg_num[2] #I wanna get rows
out_csv %>% slice(iwgr) %>% write_csv(arg[1], '_', arg[2], '_samples_exomes.csv')
printed_fams <- out_csv %>% slice(iwgr) %>% pull(patient) %>% unique()
cat(c('Записана таблица с фамилиями:\n', paste(printed_fams, collapse = '\n')))
ids_tab %>% filter(patient %in% printed_fams) %>% pull(norm_id) %>% paste0(., rep(c('_1.fq.gz', '_2.fq.gz'), each=4)) %>% write(paste0('norm_ids_', arg[1],'_', arg[2], '.txt'))
ids_tab %>% filter(patient %in% printed_fams) %>% pull(tumor_id) %>% paste0(., rep(c('_1.fq.gz', '_2.fq.gz'), each=4)) %>% write(paste0('tumor_ids_', arg[1],'_', arg[2], '.txt'))
cat(c('Получены id по фамилиям:\n', paste(printed_fams, collapse = '\n')))

