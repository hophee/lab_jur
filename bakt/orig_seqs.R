library(Biostrings)
library(tidyverse)

narna <- readDNAStringSet('all_naRNA.fasta')
narna <- narna %>% as.data.frame() %>% rownames_to_column('seq_name') %>% dplyr::rename('seq' = 'x')

red_narna <- narna %>%
  mutate(original_order = row_number()) %>%
  group_by(seq) %>%
  mutate(new_name = if(n() > 1) paste0(substring(seq, 1, 5),'_', as.character(n()),'_', as.character(nchar(seq)),'_', substring(seq, nchar(seq)-4), '-', as.character(floor(runif(1, 10000, 99999))), ' Nucleoid-associated noncoding RNA 4 (CssrE)') else seq_name) %>%
  ungroup() %>%
  arrange(original_order) %>%
  dplyr::select(-original_order)
paste0('>', red_narna$new_name, '\n', red_narna$seq) %>% write('changed_all_naRNA.fasta')
narna_uniq <- red_narna %>% 
  distinct(new_name, .keep_all = TRUE)
paste0('>', narna_uniq$new_name, '\n', narna_uniq$seq) %>% write('changed_uniq_naRNA.fasta')