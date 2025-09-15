library(tidyverse)

tablo <- readRDS("all_narna.Rds")
imp <- c('ECZV_18330','ECZV_18305','ECZV_18395','ECZV_20640','ECZV_19270')
targets <- tablo %>% filter(seq_name %in% paste(imp, 'Nucleoid-associated noncoding RNA 4 (CssrE)')) %>% pull(new_name)

bakts <- tablo %>% filter(new_name %in% targets) %>% pull(seq_name) %>% gsub(' Nucleoid-associated noncoding RNA 4 \\(CssrE\\)', '', .) %>% gsub('_\\d+$', '', .) %>%  unique()
#tablo %>% filter(new_name %in% targets) %>% pull(seq_name) %>%  gsub(' Nucleoid-associated noncoding RNA 4 \\(CssrE\\)', '', .) %>% gsub('_\\d+$', '', .) %>% table()

contain_dir <- dir('bakta_isolates')[dir('bakta_isolates') %in% paste0('bakta_annotation_',bakts)]
contained_tsv <- list.files(file.path('bakta_isolates', contain_dir), full.names = TRUE) %>% grep('\\d+\\.tsv$', ., value = T)
sapply(contained_tsv, function(x){
  read_tsv(x) %>% head()
  })
