library(tidyverse)

neighbor_of_canon <- function(nd, tab=eczv){
  idx_canon <- which((tab$gene == 'naRNA4')&(tab$is_canon==1))
  idx_non_canon <- which((tab$gene == 'naRNA4')&(tab$is_canon==0))
  can <- tab %>% dplyr::slice(c(outer((1:nd),idx_canon, '+'), outer((-nd:-1),idx_canon, '+')) %>% as.vector()) %>% mutate(canon_neighb=1)
  no_can <- tab %>% dplyr::slice(c(outer((1:nd),idx_non_canon, '+'), outer((-nd:-1),idx_non_canon, '+')) %>% as.vector()) %>% mutate(canon_neighb=0)
  neighbors <- rbind(can, no_can) %>% dplyr::arrange(start) #%>% distinct(start, .keep_all = TRUE)
  return(neighbors)
}

eczv <- read_tsv('bakt_files/ECZV.tsv', skip=2) %>% janitor::clean_names()
locus_tag <- eczv %>% filter(gene=='naRNA4') %>% pull(locus_tag) %>% sort()
is_canon <- c(0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1)
canon_tab <- data.frame(locus_tag, is_canon) %>% slice(-58) #в 58 строке спорный образец
all_neighb <- left_join(eczv, canon_tab, by='locus_tag') %>% neighbor_of_canon(1, tab=.) %>% select(start, type,canon_neighb, gene, locus_tag, is_canon)

non_can <- all_neighb %>% filter(type %in% c('cds', 'ncRNA')) %>% filter(canon_neighb==0) %>% pull(type) %>% as.factor()
can <- all_neighb %>% filter(type %in% c('cds', 'ncRNA')) %>% filter(canon_neighb==1) %>% pull(type) %>% as.factor()
contingency_table <- table(group = c(rep("non_can", length(non_can)), rep("can", length(can))),
                           value = c(non_can, can))
chi2_test <- chisq.test(contingency_table)
#Значение p-value
print(chi2_test$p.value)
#Наблюдаемые значения
print(chi2_test$observed)