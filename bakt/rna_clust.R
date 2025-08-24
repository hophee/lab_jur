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

for_vis_df <- bind_rows(cl_sz_plus_strand, cl_sz_minus_strand) %>% 
  select(-grep('*\\.y', colnames(cl_sz_minus_strand))) %>% dplyr::rename('end'='stop.x')
colnames(for_vis_df) <- sub('\\.x', '', colnames(for_vis_df))
narna_vis <- for_vis_df %>% filter(gene=='naRNA4') %>% select(1, 3, 4, 5, 17)

plus_vis <- cl_sz_plus_strand %>% select(-grep('*\\.y', colnames(cl_sz_plus_strand))) %>% dplyr::rename('end'='stop.x') %>% filter(gene.x=='naRNA4') %>%  select(1, 3, 4, 5, 17)
colnames(plus_vis) <- sub('\\.x', '', colnames(plus_vis))

min_vis <- cl_sz_minus_strand %>% select(-grep('*\\.y', colnames(cl_sz_minus_strand))) %>% dplyr::rename('end'='stop.x') %>% filter(gene.x=='naRNA4') %>%  select(1, 3, 4, 5, 17)
colnames(min_vis) <- sub('\\.x', '', colnames(min_vis))

chromosome.size <- max(eczv$stop)
circos.par(cell.padding = c(0, 0, 0, 0))  # Убрать отступы
circos.initialize("contig_1", xlim = c(0, chromosome.size))

strand_colors <- c("+" = "darkgreen", "-" = "darkred")

circos.track(
  ylim = c(0, 1),
  track.height = 0.3,  # Высота дорожки
  panel.fun = function(x, y) {
    if (nrow(min_vis) == 0) return()
    
    for (i in 1:nrow(min_vis)) {
      start <- min_vis$start[i]
      end <- min_vis$end[i]
      cluster_size <- min_vis$cluster_size[i]
      
      # Рисуем прямоугольник
      circos.rect(start, 0, end, 1, 
                  col = strand_colors["-"], border = strand_colors["-"], lwd = 1)
      
      
        center_x <- (start + end) / 2
        circos.text(center_x, 0.5, labels = as.character(cluster_size),
                    facing = "bending.inside", niceFacing = TRUE,
                    cex = 0.7, col = "black", font = 2)
    }
  }
)

circos.track(
  ylim = c(0, 1),
  track.height = 0.3,  # Высота дорожки
  panel.fun = function(x, y) {
    if (nrow(plus_vis) == 0) return()
    
    for (i in 1:nrow(plus_vis)) {
      start <- plus_vis$start[i]
      end <- plus_vis$end[i]
      cluster_size <- plus_vis$cluster_size[i]
      
      # Рисуем прямоугольник
      circos.rect(start, 0, end, 1, 
                  col = strand_colors["+"], border = strand_colors["+"], lwd = 0.5)
      
      # Добавляем надпись (логика такая же как раньше)
      show_label <- TRUE
      if (cluster_size > 1) {
        # Проверяем, является ли это первым в кластере
        if (i > 1 && plus_vis$start[i] <= plus_vis$end[i-1]) {
          show_label <- FALSE
        }
      }
      
      if (show_label) {
        center_x <- (start + end) / 2
        circos.text(center_x, 0.5, labels = as.character(cluster_size),
                    facing = "bending.inside", niceFacing = TRUE,
                    cex = 0.7, col = "black", font = 2)
      }
    }
  }
)
circos.clear()


