library(tidyverse)
library(circlize)

eczv <- read_tsv('bakt_files/ECZV.tsv', skip=2) %>% janitor::clean_names()
dif_exp <- readxl::read_xlsx('rna_dif_exp.xlsx', sheet = 1) %>% janitor::clean_names() %>% dplyr::rename('locus_tag' = 'gene_2')

full_tab <- left_join(eczv, dif_exp, by='locus_tag')
plus_strand <- full_tab %>% filter(strand.x == '+')
minus_strand <- full_tab %>% filter(strand.x == '-')


calculate_cluster_size <- function(strand_df) {
  strand_df %>%
    group_by(number_sequence_id) %>%
    arrange(
      if (first(strand_df$strand.x) == "+") start.x else -start.x
    ) %>%
    mutate(
      cluster_size = {
        rle_res <- rle(gene.x)
        rep(rle_res$lengths, rle_res$lengths)
      }
    ) %>%
    ungroup()
}
cl_sz_plus_strand <- calculate_cluster_size(plus_strand)
cl_sz_minus_strand <- calculate_cluster_size(minus_strand)

rle1 <- cl_sz_plus_strand %>% filter(gene.x=='naRNA4') %>% pull(cluster_size) %>% rle()
rep(rle1$values, rle1$lengths/rle1$values) %>% 
  qplot(geom = "bar", fill=I("lightblue"), color = I("black"))+
  theme_minimal()+
  scale_y_continuous(breaks = seq(1, max(table(rep(rle1$values, rle1$lengths/rle1$values))), by = 2)) +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())+
  ggtitle('Cluster size distribution on the plus strand')+
  xlab('Cluster size')
rle2 <- cl_sz_minus_strand %>% filter(gene.x=='naRNA4') %>% pull(cluster_size) %>% rle()
rep(rle2$values, rle2$lengths/rle2$values) %>% 
  qplot(geom = "bar", fill=I("lightblue"), color = I("black"))+
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())+
  scale_y_continuous(breaks = seq(1, max(table(rep(rle2$values, rle2$lengths/rle2$values))), by = 2)) +
  ggtitle('Cluster size distribution on the minus strand')+
  xlab('Cluster size')

for_vis_df <- bind_rows(cl_sz_plus_strand, cl_sz_minus_strand) %>% select(-grep('*\\.y', colnames(cl_sz_minus_strand)))
colnames(for_vis_df) <- sub('\\.x', '', colnames(cl_sz_minus_strand))

chromosome.size <- max(eczv$stop)
circos.par(cell.padding = c(0, 0, 0, 0))  # Убрать отступы
circos.initialize("chr1", xlim = c(0, chromosome.size))

