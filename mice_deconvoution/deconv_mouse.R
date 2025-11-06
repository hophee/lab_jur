library(celldex)
library(tidyverse)
library(SingleCellExperiment)
library(FARDEEP)

ref <- fetchReference('mouse_rnaseq', fetchLatestVersion('mouse_rnaseq'))
mouse_atlas <- assay(ref, 'logcounts')
genes_atlas <- rownames(mouse_atlas)
atlas_tpm <- 2^mouse_atlas - 1
rownames(atlas_tpm) <- genes_atlas
cell_types <- colData(ref)$label.main


counts <- read_tsv("salmon.merged.gene_tpm.tsv")
counts <- counts %>% group_by(gene_name) %>% summarise(mice_11_control_PA = sum(mice_11_control_PA), mice_17_exp_Ecoli_PA = sum(mice_17_exp_Ecoli_PA)) %>% column_to_rownames(var='gene_name')

comon_genes <- intersect(rownames(atlas_tpm), rownames(counts))
atlas_tpm <- atlas_tpm[comon_genes,]
counts <- counts[comon_genes,]

res <- fardeep(X = atlas_tpm, Y = counts, lognorm = TRUE, QN = FALSE)
