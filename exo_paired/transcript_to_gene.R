library(biomaRt)
library(tidyverse)

abundance <- read_tsv("abundance.tsv") %>% filter(tpm != 0)
#transcript_ids <- strsplit(abundance$target_id, split = '|', fixed = T) %>% sapply('[', 1)
abundance$target_id <- gsub('\\..+$', '', abundance$target_id)

mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
symbol <- getBM(filters = "ensembl_transcript_id",
                attributes = c("ensembl_transcript_id","external_gene_name"),
                values = abundance$target_id, 
                mart = mart)
abundance <- abundance %>% mutate(target_id = symbol)
new <- merge(abundance, symbol, by.x='target_id', by.y='ensembl_transcript_id') %>% dplyr::select(external_gene_name, tpm) 
colnames(new) <- c('TargetID', 'TPM')
write_tsv(new, file='annotated_rna_expression.tsv')
