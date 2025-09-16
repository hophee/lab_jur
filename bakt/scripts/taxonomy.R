library(tidyverse)
library(rentrez)
library(taxize)
library(ggtree)
library(ggthemes)

taxonomy <- read_csv('9534URHV013-Alignment-HitTable.csv', 
                     col_names=c('query acc.ver', 'subject acc.ver', '% identity', 'alignment length', 'mismatches', 'gap opens', 'q. start', 'q. end', 's. start', 's. end', 'evalue', 'bit score')) %>% 
  janitor::clean_names()
uni_ids <- unique(taxonomy$subject_acc_ver)

get_acc <- sapply(uni_ids, function(x) {
  res <- entrez_search(db='nucleotide', term = paste0(x, "[Accession]"), retmax=1)
  if(length(res$ids) == 0) return(NA)
  return(res$ids[1])
})
unioned_data <- data.frame(uni_ids, unlist(get_acc)) %>% rename('GI'='unlist.get_acc.', 'subject_acc_ver' = 'uni_ids')

step <- 42
res_names <- vector()
res_taxid <- vector()
for (i in seq(1, length(unioned_data$GI), step)) {
  end <- min(i+step-1, length(unioned_data$GI))
  get_names <- sapply(unioned_data$GI[i:end], function(x){entrez_summary(db = "nucleotide", id = x)$organism})
  get_taxid <- sapply(unioned_data$GI[i:end], function(x){entrez_summary(db = "nucleotide", id = x)$taxid})
  res_names <- c(res_names, get_names)
  res_taxid <- c(res_taxid, get_taxid)
  Sys.sleep(0.2)
}
unioned_data <- unioned_data %>% mutate('org_name' = res_names, 'taxid'=res_taxid)
taxonomy <- left_join(taxonomy, unioned_data, by='subject_acc_ver')
taxonomy <- taxonomy %>% mutate('cleaned_names'= sapply(strsplit(taxonomy$org_name, ' '), function(x){ paste(x[1], x[2])}))
#strange
clas <- classification(unioned_data$taxid, db='ncbi')
tre <- class2tree(unique(clas))
ggtree(tre$phylo)+
  geom_tiplab()+
  theme_tree2()

taxonomy_counts <- taxonomy %>%
  group_by(query_acc_ver, org_name) %>%
  summarise(count = n(), .groups = "drop")

taxonomy_counts_norm <- taxonomy %>% filter(evalue<10^(-5)) %>% 
  group_by(query_acc_ver, cleaned_names) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(query_acc_ver) %>%
  mutate(freq = count / sum(count)) %>%
  ungroup() %>% filter(count > 5,freq>0.1)

tax_filtered_summed <- taxonomy_counts_norm %>% group_by(cleaned_names) %>% summarise(all_counts=sum(count), median_freq=median(freq), .groups = "drop") %>% 
  arrange(desc(all_counts)) %>%   mutate(cleaned_names = factor(cleaned_names, levels = cleaned_names[order(-all_counts)])) %>% filter(cleaned_names!='uncultured bacterium')

pdf('spicies.pdf', width=10, height = 8, family = "Times")
ggplot(data=tax_filtered_summed)+
  geom_point(aes(y=log10(all_counts), x=median_freq, color=cleaned_names, size=10))+
  geom_text(aes(y = log10(all_counts), x = median_freq, label = cleaned_names), vjust = -0.5, size = 3)+
  scale_color_brewer(palette = "Spectral")+
  xlab('Медианная частота встерчаемости')+
  labs(y='Десятичный логарифм числа вхождений')+
  theme(legend.position = "none")
ggplot()+
  geom_col(data=tax_filtered_summed, aes(y=all_counts, x=cleaned_names, fill=cleaned_names))+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position = "none")+
  labs(y='Суммарное число вхождений', x='')+
  scale_fill_brewer(palette = "Spectral")
ggplot()+
  geom_col(data=tax_filtered_summed, aes(y=median_freq, x=cleaned_names, fill=cleaned_names))+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")+
  labs(y='Медианная частота встречаемости', x='')+
  scale_fill_brewer(palette = "Spectral")
dev.off()


tree <- read.tree('cleaned_nucleo_narna.phb')
ggtree(tree) + geom_tiplab() + theme_tree2()
