library(tidyverse)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)

data <- read_delim('tumor_with_rg.call.cns') %>% filter(abs(log2) > 0.5)

gr <- GRanges(seqnames=data$chromosome, ranges=IRanges(start=data$start, end=data$end), log2=data$log2)


genes <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene)
hits <- findOverlaps(gr, genes)
annotated <- cbind(data[queryHits(hits),], 
                   gene_id=names(genes[subjectHits(hits)]))
annotated$gene <- mapIds(org.Hs.eg.db,
                           keys=annotated$gene_id,
                           column="SYMBOL",
                           keytype="ENTREZID")

annotated %>% pull(gene) %>% unique() %>% write(file='high_cnv_gene_list.txt')
