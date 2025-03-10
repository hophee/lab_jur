library(tidyverse)
library(janitor)
library(gtools)
library(purrr)

#оценка распределения качества неотфильтрованных вариантов
qual <- scan('qual_values.txt', numeric())
ggplot() +
  geom_histogram(aes(x=qual), bins = 1000)

#функция для грепа по множеству паттернов
find_patterns <- function(pattern_vector, df, column_name) {
  if (!column_name %in% colnames(df)) {
    stop("Указанная колонка не найдена в датафрейме")
  }
  
  result <- map_lgl(df[[column_name]], ~ any(map_lgl(pattern_vector, ~ grepl(.x, .y, ignore.case = T), .y = .x)))
  
  return(result)
}

anovar <- read_csv('annovar.csv', na = c("", "NA", '.')) %>% clean_names() %>% mutate(chr = factor(chr)) %>%  mutate(chr = factor(chr, levels = unique(chr)[mixedorder(chr)])) %>%  
  select(where(~ !all(is.na(.))))
vep <- read_delim('vep.txt', na = c("", "NA", '.')) %>% clean_names() %>%  select(where(~ !all(is.na(.)))) %>% mutate(gnom_a_de_af = as.numeric(gnom_a_de_af))

#PVS1 варианты
anovar_pvs1 <- anovar %>% filter(find_patterns(c('frameshift', 'stopgain', 'stoploss', 'startloss'), anovar, 'exonic_func_ref_gene_with_ver')) %>% filter(gnomad41_exome_af < 0.01)
vep_pvs1 <- vep %>% filter(find_patterns(c('frameshift_variant', 'stop_gained', 'stop_lost', 'start_lost', 'splice_acceptor_variant', 'splice_donor_variant', 'transcript_ablation'), vep, 'consequence')) %>% 
  filter(gnom_a_de_af < 0.001) %>% filter(!grepl('benign', clin_sig, ignore.case = T)) %>% filter(!grepl('benign', poly_phen, ignore.case = T))
vep_genes <- vep_pvs1 %>% pull(symbol) %>% unique()

pvs1_vep_ano_intersect <- anovar_pvs1 %>% filter(gene_ref_gene_with_ver %in% vep_genes) 

#SNV варианты
most_pat_anovar_pp3 <- anovar %>% filter(exonic_func_ref_gene_with_ver %in% 'nonsynonymous SNV') %>% filter(sift_score < 0.05) %>% 
  filter(polyphen2_hdiv_score > 0.9) %>%  filter(provean_score < -2.5) %>% filter(meta_svm_score > 0) %>% filter(cadd_phred > 20) %>% 
  filter(gnomad41_exome_af < 0.001) %>% filter(!(clnsig %in% c('Likely_benign', 'Benign', 'Benign/Likely_benign')))
vep_path_snps <- vep %>% filter(consequence %in% 'missense_variant') %>% filter(gnom_a_de_af < 0.001) %>% 
  filter(!grepl('benign', poly_phen, ignore.case=T)) %>% filter(!grepl('tolerated', sift, ignore.case=T)) %>% 
  filter(!grepl('benign', clin_sig, ignore.case=T)) %>% filter(as.numeric(cadd_phred)>30)

vep_snips_genes <- vep_path_snps %>% pull(symbol) %>% unique()
pp3_intersect <- most_pat_anovar_pp3 %>% filter(gene_ref_gene_with_ver %in% vep_snips_genes) 


#Полученные варианты должны быть в количестве, которое реалистично рассмотреть вручную
#Если вариантов много -> фильтровать жёстче
#Если вариантов мало/в них ничего не найдено -> изменить параметры в сторону менее жёсткой фильтрации
most_pathogenic <- rbind(pp3_intersect, pvs1_vep_ano_intersect) %>% select(where(~ !all(is.na(.))))