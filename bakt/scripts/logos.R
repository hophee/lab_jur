library(ggseqlogo)
library(ggtree)
library(tidyverse)
library(reshape2)
library(ape)


make_logo <- function(mult_align, name, st='dna', mt='bits'){
  aligned_file <- readDNAMultipleAlignment(mult_align)
  aligned_mat <- as.matrix(aligned_file)
  cons_mat <- consensusMatrix(aligned_mat)
  ggplot()+
    geom_logo(data=cons_mat, seq_type = st, method = mt)+
    theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5))+
    scale_x_continuous(breaks = c(seq(0, dim(cons_mat)[2], 5), dim(cons_mat)[2]))+
    labs(title=name)
}

make_logo("all_naRNA.aln", 'Выравниванивание всех последовательностей')
make_logo("czv_naRNA.aln", 'Выравниванивание по ECZV')
make_logo('changed_uniq_naRNA.aln', 'Уникальные последовательности')

#??????
dnd_text <- readLines("all_naRNA.dnd")
dnd_text <- paste(dnd_text, collapse = "")
dnd_text <- gsub("\\s+", "", dnd_text)
#dnd_text <- gsub(":-?\\d+\\.\\d+", ":0.00001", dnd_text)
tree <- read.tree(text = dnd_text)
ggtree(tree) + 
  geom_tiplab(size = 0.5) + 
  ggtitle("Phylogenetic Tree from .dnd")+
  coord_flip()
