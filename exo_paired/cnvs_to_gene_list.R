library(tidyverse)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)

data <- read_delim('tumor_sample_with_rg.somatic.call.cns') %>% filter(chromosome != 'chrY') %>%  mutate(width = end-start, effect = width* abs(log2)) 

gr <- GRanges(seqnames=data$chromosome, ranges=IRanges(start=data$start, end=data$end), log2=data$log2)


genes <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene)
hits <- findOverlaps(gr, genes)
annotated <- cbind(data[queryHits(hits),], 
                   gene_id=names(genes[subjectHits(hits)]))
annotated$gene <- mapIds(org.Hs.eg.db,
                           keys=annotated$gene_id,
                           column="SYMBOL",
                           keytype="ENTREZID")


dels <- annotated %>% filter(log2< (-0.8))
deleted_genes <- dels %>% pull(gene) %>% unique()
deleted_genes[grep('PTP', deleted_genes)]

ampli <- annotated %>% filter(log2>0.585)
ampli <- ampli[-grep('LOC', ampli$gene),]
ampli <- ampli[-grep('LINC', ampli$gene),]
ampli_genes <- ampli %>% pull(gene) %>% unique()
ampli_genes[grep('RB', ampli_genes)]
write(ampli_genes, file='amplified_genes.txt')

#same
annotated[grep('CD274', annotated$gene),]
annotated %>% filter(gene=='CD274')

hard_filtered <- annotated %>% filter(effect > quantile(unique(effect), .95))

annotated %>% pull(gene) %>% unique() %>% sort() %>% write(file='very_high_cnv_gene_list.txt')

ggplot() +
  geom_point(data=data, aes(x=effect, y=weight))
annotated %>% filter(chromosome == 'chrX') %>% pull(gene) %>% write(file='chrX_list.txt')
