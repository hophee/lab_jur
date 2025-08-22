library(tidyverse)

eczv <- read_tsv('bakt_files/ECZV.tsv', skip=2) %>% janitor::clean_names()
dif_exp <- readxl::read_xlsx('rna_dif_exp.xlsx', sheet = 1) %>% janitor::clean_names() %>% dplyr::rename('locus_tag' = 'gene_2')

full_tab <- left_join(eczv, dif_exp, by='locus_tag')
plus_strand <- full_tab %>% filter(strand.x == '+')
minus_strand <- full_tab %>% filter(strand.x == '-')

which(plus_strand$gene.x=='naRNA4')

plus_strand %>% 
  mutate(
    cluster_id = cumsum(gene.x != lag(gene.x, default = first(gene.x)))
  ) %>%
  group_by(cluster_id, .add = TRUE) %>%
  mutate(
    cluster_size = n()
  ) %>%
  ungroup() %>%
  select(-cluster_id) %>% filter(gene.x=='naRNA4')

calculate_cluster_size <- function(strand_df) {
  strand_df %>%
    group_by(number_sequence_id) %>%
    arrange(
      if (first(strand_df$strand.x) == "+") start.x else -start.x
    ) %>%
    # Используем rle для надежного определения кластеров
    mutate(
      cluster_size = {
        # Преобразуем в character, чтобы избежать проблем с факторами
        gene_vec <- as.character(gene.x)
        rle_res <- rle(gene_vec)
        rep(rle_res$lengths, rle_res$lengths)
      }
    ) %>%
    ungroup()
}

calculate_cluster_size(plus_strand) %>% filter(gene.x=='naRNA4')
