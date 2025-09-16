library(ggbio)
library(GenomicRanges)
library(dplyr)

# Преобразуем данные в GRanges
narna_gr <- makeGRangesFromDataFrame(
  narna_vis,
  seqnames.field = "number_sequence_id",
  start.field = "start",
  end.field = "end",
  strand.field = "strand"
)


# Добавляем cluster_size как метаданные
mcols(narna_gr)$cluster_size <- narna_vis$cluster_size

genome_circle <- GRanges(
  seqnames = "contig_1",  # или ваше имя последовательности
  IRanges(start = 1, end = max(eczv$stop)+100),
  strand = "*"
)
# Создаем отдельные объекты для каждого strand
plus_gr <- narna_gr[strand(narna_gr) == "+"]
minus_gr <- narna_gr[strand(narna_gr) == "-"]


sig_gr <- full_tab %>% filter(significant=='yes') %>% filter(!is.na(gene.x))%>% 
  makeGRangesFromDataFrame(. ,
                           seqnames.field = "gene.x",
                           start.field = "start.x",
                           end.field = "stop.x",
                           strand.field = "strand.x")


# Создаем график с разными цветами для каждого strand
p <- ggbio()+
  circle(plus_gr, geom='rect', aes(fill=strand, col=strand))+
  circle(minus_gr, geom='rect', aes(fill=strand, col=strand))+
  circle(plus_gr, geom='text', aes(label=cluster_size, vjust = 12))+
  circle(minus_gr, geom='text', aes(label=cluster_size, vjust = -0.25))+
  
  circle(genome_circle, geom = "rect", 
         aes(fill = "background"), 
         fill = "lightgray", alpha = 0.01, size = 0.5)+
  labs(fill='Strand of naRNA4', title='Расположение naRNA4 по геному')+
  theme(plot.title = element_text(hjust = 0.5))+
  guides(col=F)


# Добавляем надписи отдельно
print(p)



#-----------------------
bin_size <- 10000
chrom_length <- max(eczv$stop)
bins <- tileGenome(seqlengths = c(contig_1 = chrom_length), tilewidth = bin_size, cut.last.tile.in.chrom = TRUE)
mcols(bins)$density <- countOverlaps(bins, narna_gr)
genome_circle <- GRanges(seqnames = "contig_1", IRanges(start = 1, end = chrom_length))
p_density <- ggbio() +
  circle(genome_circle, geom = "rect",
         fill = "lightgray", alpha = 0.1, size = 0.5) +
  circle(bins, geom = "rect", aes(fill = density), color = NA) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  labs(title = "Карта плотности naRNA4 по геному", fill='Число генов\nв окне 10кб') +
  theme(plot.title = element_text(hjust = 0.5))

print(p_density)

#-----------

min_dist <- function(x1, x2, y1, y2){
  dists1 <- x1 - y1
  dists2 <- x2 - y2
  closest_idx1 <- which.min(abs(dists1))
  closest_idx2 <- which.min(abs(dists2))
  if (abs(dists1[closest_idx1])<abs(dists2[closest_idx2])){
    return(dists1[closest_idx1])
  }else{
    return(dists2[closest_idx2])
  }
}
min_dist2 <- function(x1, x2, y1, y2, df1, df2){
  dists1 <- x1 - y1
  dists2 <- x2 - y2
  closest_idx1 <- which.min(abs(dists1))
  closest_idx2 <- which.min(abs(dists2))
  if (abs(dists1[closest_idx1])<abs(dists2[closest_idx2])){
    min_res <- dists1[closest_idx1]
    idx <- closest_idx1
    }else{
    min_res <- dists2[closest_idx2]
    idx <- closest_idx2
    }
  lct <- df1$locus_tag[closest_idx2]
  cl_sz <- df2$cluster_size[which(df2$locus_tag==lct)]
  return(c(min_res, cl_sz))
}
narna_plus <- full_tab %>% filter((gene.x=='naRNA4') & (strand.x=='+'))
narna_minus <- full_tab %>% filter((gene.x=='naRNA4') & (strand.x=='-'))
sig_plus <- full_tab %>% filter((strand.x=='+') & (significant=='yes')) %>% 
  mutate(rna_dist=mapply(function(x1, x2){min_dist2(x1, x2, narna_plus$stop.x, narna_plus$start.x, narna_plus, cl_sz_plus_strand)[1]}, start.x, stop.x)) %>% 
  mutate(cl_sz=mapply(function(x1, x2){min_dist2(x1, x2, narna_plus$stop.x, narna_plus$start.x, narna_plus, cl_sz_plus_strand)[2]}, start.x, stop.x))
sig_minus <- full_tab %>% filter((strand.x=='-') & (significant=='yes')) %>% 
  mutate(rna_dist=mapply(function(x1, x2){min_dist2(x1, x2, narna_minus$stop.x, narna_minus$start.x, narna_minus, cl_sz_minus_strand)[1]}, start.x, stop.x)) %>% 
  mutate(cl_sz=mapply(function(x1, x2){min_dist2(x1, x2, narna_minus$stop.x, narna_minus$start.x, narna_minus, cl_sz_minus_strand)[2]}, start.x, stop.x))


plot_data <- rbind(sig_plus, sig_minus) %>% mutate(log2fold_change=as.numeric(log2fold_change))
ggplot()+
  geom_point(data=plot_data, aes(y=log2fold_change, x=log10(abs(rna_dist)), col=cl_sz))+
  scale_color_gradient(low='blue', high='red')+
  labs(title = 'Отношение расстояния от ближайшего naRNA4 к log2(FC)', color='Размер \nближайшего \nкластера', 
       x='log10(|Расстояние до naRNA4|)', y='log2(FC)')+
  theme_minimal()+
  theme(plot.title = element_text(hjust=0.5))+
  ggrepel::geom_label_repel(data=(plot_data %>% filter(log10(abs(rna_dist))<3.3)), 
             aes(label=gene.x, y=log2fold_change+1.5, x=log10(abs(rna_dist)), alpha=0.05))+
  guides(alpha=F)

mod <- plot_data %>% filter(log10(abs(rna_dist))<3.3) %>% filter(log2fold_change>2)
ggplot()+
  geom_point(data=mod, aes(x=log2fold_change, y=log10(abs(rna_dist)), col=cl_sz))+
  scale_color_gradient(low='blue', high='red')+
  labs(title = 'Отношение расстояния от ближайшего naRNA4 к log2(FC)', color='Размер \nближайшего \nкластера', 
       y='log10(|Расстояние до naRNA4|)', x='log2(FC)')+
  theme_minimal()+
  theme(plot.title = element_text(hjust=0.5))+
  ggrepel::geom_label_repel(data=mod, 
                            aes(label=gene.x, x=log2fold_change, y=log10(abs(rna_dist)), alpha=0.05))+
  guides(alpha=F)+
  xlim(c(0, 15))+
  ylim(c(0, 3.5))+
  geom_vline(xintercept =2, color='red', linetype=2)
