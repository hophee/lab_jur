library(tidyverse)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)

data <- read_delim('tumor_with_rg.call.cns') %>% mutate(width = end-start, effect = width* abs(log2)) %>% filter(chromosome != 'chrY')

gr <- GRanges(seqnames=data$chromosome, ranges=IRanges(start=data$start, end=data$end), log2=data$log2)


genes <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene)
hits <- findOverlaps(gr, genes)
annotated <- cbind(data[queryHits(hits),], 
                   gene_id=names(genes[subjectHits(hits)]))
annotated$gene <- mapIds(org.Hs.eg.db,
                           keys=annotated$gene_id,
                           column="SYMBOL",
                           keytype="ENTREZID")

hard_filtered <- annotated %>% filter(effect > quantile(unique(effect), .95))
annotated %>% pull(gene) %>% unique() %>% sort() %>% write(file='very_high_cnv_gene_list.txt')

ggplot() +
  geom_point(data=data, aes(x=effect, y=weight))
annotated %>% filter(chromosome == 'chrX') %>% pull(gene) %>% write(file='chrX_list.txt')
