library(tidyverse)
library(clusterProfiler)
library(DOSE)
library(enrichplot)
library(org.EcK12.eg.db)

eczv <- read_tsv('bakt_files/ECZV.tsv', skip=2) %>% janitor::clean_names()
ids <- read_lines('narna_ids.txt')

eczv <- eczv %>% 
  mutate(db_xrefs = strsplit(db_xrefs, ",\\s*")) %>% 
  unnest(db_xrefs) %>%
  separate(db_xrefs, into = c("db", "id"), sep = ":", extra = "merge") %>% 
  distinct(number_sequence_id, type, start, stop, strand, locus_tag, gene, product, db, .keep_all = TRUE) %>% 
  pivot_wider(
    names_from = db, 
    values_from = id,
    values_fill = NA,
    values_fn = ~ paste(unique(.x), collapse = "; ")  # Объединяем дубли через "; "
  )

#nd - neighbor distance
get_plots <- function(nd, pcut=1, qcut=1){
  idxes <- which(eczv$gene == 'naRNA4')
  right_neighbors <- eczv %>% dplyr::slice(outer((1:nd),idxes, '+') %>% as.vector()) %>% dplyr::filter(gene!='naRNA4') %>% mutate(side='right')
  left_neighbors <- eczv %>% dplyr::slice(outer((-nd:-1),idxes, '+') %>% as.vector()) %>% dplyr::filter(gene!='naRNA4')%>% mutate(side='left')
  neighbors <- rbind(right_neighbors, left_neighbors) %>% dplyr::arrange(start)
  go_right <- neighbors %>% dplyr::filter(side=='right') %>% pull(GO) %>% na.omit() %>% paste0('GO:', .)
  go_left <- neighbors %>% dplyr::filter(side=='left') %>% pull(GO) %>% na.omit() %>% paste0('GO:', .)
  ego_right <- enrichGO(gene = go_right, 
                        OrgDb = org.EcK12.eg.db,
                        keyType = "GO", 
                        ont = "ALL",
                        pvalueCutoff = pcut,
                        pAdjustMethod = "BH",
                        qvalueCutoff = qcut,
                        readable = TRUE)
  right_dots <- dotplot(ego_right, showCategory = 10, 
          font.size = 10, 
          color = "p.adjust", 
          split = "ONTOLOGY", 
          title = paste0('Right genes, nd=', as.character(nd))) + 
    facet_grid(ONTOLOGY~., scales = "free")
  right_go_bp <- enrichGO(gene = go_right, 
           OrgDb = org.EcK12.eg.db,
           keyType = "GO", 
           ont = "BP",
           pvalueCutoff = pcut,
           pAdjustMethod = "BH",
           qvalueCutoff = qcut,
           readable = TRUE) %>% goplot() +
    ggtitle(paste0('Right genes, Biological Process onthology, nd=', as.character(nd)))+
    theme(plot.title = element_text(hjust = 0.5))
  right_go_cc <- enrichGO(gene = go_right, 
                          OrgDb = org.EcK12.eg.db,
                          keyType = "GO", 
                          ont = "CC",
                          pvalueCutoff = pcut,
                          pAdjustMethod = "BH",
                          qvalueCutoff = qcut,
                          readable = TRUE) %>% goplot() +
    ggtitle(paste0('Right genes, Cellular Component onthology, nd=', as.character(nd)))+
    theme(plot.title = element_text(hjust = 0.5))
  right_go_MF <- enrichGO(gene = go_right, 
                          OrgDb = org.EcK12.eg.db,
                          keyType = "GO", 
                          ont = "MF",
                          pvalueCutoff = pcut,
                          pAdjustMethod = "BH",
                          qvalueCutoff = qcut,
                          readable = TRUE) %>% goplot() +
    ggtitle(paste0('Right genes, Molecular Function onthology, nd=', as.character(nd)))+
    theme(plot.title = element_text(hjust = 0.5))
  
  ego_left <- enrichGO(gene = go_left, 
                        OrgDb = org.EcK12.eg.db,
                        keyType = "GO", 
                        ont = "ALL",
                        pvalueCutoff = pcut,
                        pAdjustMethod = "BH",
                        qvalueCutoff = qcut,
                        readable = TRUE)
  left_dots <- dotplot(ego_left, showCategory = 10, 
                        font.size = 10, 
                        color = "p.adjust", 
                        split = "ONTOLOGY", 
                        title = paste0('Left genes, nd=', as.character(nd))) + 
    facet_grid(ONTOLOGY~., scales = "free")
  left_go_bp <- enrichGO(gene = go_left, 
                          OrgDb = org.EcK12.eg.db,
                          keyType = "GO", 
                          ont = "BP",
                          pvalueCutoff = pcut,
                          pAdjustMethod = "BH",
                          qvalueCutoff = qcut,
                          readable = TRUE) %>% goplot() +
    ggtitle(paste0('Left genes, Biological Process onthology, nd=', as.character(nd)))+
    theme(plot.title = element_text(hjust = 0.5))
  left_go_cc <- enrichGO(gene = go_left, 
                          OrgDb = org.EcK12.eg.db,
                          keyType = "GO", 
                          ont = "CC",
                          pvalueCutoff = pcut,
                          pAdjustMethod = "BH",
                          qvalueCutoff = qcut,
                          readable = TRUE) %>% goplot() +
    ggtitle(paste0('Left genes, Cellular Component onthology, nd=', as.character(nd)))+
    theme(plot.title = element_text(hjust = 0.5))
  left_go_MF <- enrichGO(gene = go_left, 
                          OrgDb = org.EcK12.eg.db,
                          keyType = "GO", 
                          ont = "MF",
                          pvalueCutoff = pcut,
                          pAdjustMethod = "BH",
                          qvalueCutoff = qcut,
                          readable = TRUE) %>% goplot() +
    ggtitle(paste0('Left genes, Molecular Function onthology, nd=', as.character(nd)))+
    theme(plot.title = element_text(hjust = 0.5))
  return(list(right_dots, right_go_MF, right_go_cc, right_go_bp, left_dots, left_go_MF, left_go_cc, left_go_bp))
}


get_plots(nd=1, pcut=0.05)
get_plots(nd=1, pcut=0.05)







