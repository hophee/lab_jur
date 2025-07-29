library(tidyverse)

meta <- read_tsv('metadata.tsv') %>% janitor::clean_names()
meta <- meta %>% select(1, 2, 5, 6, 7, 8) %>% group_by(familia) 
get_id_of_surname <- function(fam, nk_num){
  nk <- c("Соматическая РНК","Соматическая ДНК","Герминальная ДНК")
  ids <- meta %>% filter(familia %in% fam) %>% filter(nukleinovaa_kislota==nk[nk_num]) %>% pull(1)
  return(c(paste0(ids,'_1.fq.gz'), paste0(ids,'_2.fq.gz')))
}
meta %>% pull(2) %>% unique() %>% tail(n=4) %>% get_id_of_surname(.,nk_num=2) %>% write(file='somatic_tn4.txt')
meta %>% pull(2) %>% unique() %>% tail(n=4) %>% get_id_of_surname(.,nk_num=3) %>% write(file='germ_tn4.txt')

f4 <- meta %>% pull(2) %>% unique() %>% tail(n=4)
meta %>% filter(familia %in% f4) %>% filter(nukleinovaa_kislota %in% c("Соматическая ДНК","Герминальная ДНК")) %>% select(1, 2, 4)
