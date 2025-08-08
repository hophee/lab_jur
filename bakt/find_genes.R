library(tidyverse)
library(clusterProfiler)
library(DOSE)
library(enrichplot)
library(org.EcK12.eg.db)

eczv <- read_tsv('bakt_files/ECZV.tsv', skip=2) %>% janitor::clean_names()

eczv <- eczv %>% 
  mutate(db_xrefs = strsplit(db_xrefs, ",\\s*")) %>% 
  unnest(db_xrefs) %>%
  separate(db_xrefs, into = c("db", "id"), sep = ":", extra = "merge") %>% 
  distinct(number_sequence_id, type, start, stop, strand, locus_tag, gene, product, db, .keep_all = TRUE) %>% 
  pivot_wider(
    names_from = db, 
    values_from = id,
    values_fill = NA,
    values_fn = ~ paste(unique(.x), collapse = "; "))
entr <- bitr(
  eczv$gene,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.EcK12.eg.db) %>% rename('SYMBOL' = 'gene')
eczv <- left_join(eczv, entr, by='gene')
keggs <- bitr_kegg(eczv$ENTREZID, fromType='ncbi-geneid', toType='kegg', organism='eco', drop = TRUE) %>% rename('ncbi-geneid' = 'ENTREZID')
eczv <- left_join(eczv, keggs, by='ENTREZID')

get_neighbor <- function(nd){
  idxes <- which(eczv$gene == 'naRNA4')
  right_neighbors <- eczv %>% dplyr::slice(outer((1:nd),idxes, '+') %>% as.vector()) %>% dplyr::filter(gene!='naRNA4') %>% mutate(side='right')
  left_neighbors <- eczv %>% dplyr::slice(outer((-nd:-1),idxes, '+') %>% as.vector()) %>% dplyr::filter(gene!='naRNA4')%>% mutate(side='left')
  neighbors <- rbind(right_neighbors, left_neighbors) %>% dplyr::arrange(start) %>% distinct(start, .keep_all = TRUE)
  return(neighbors)
}

get_goplot <- function(genelist, ont_type, side_num, nd, pcut_gp=0.05, qcut_gp=1, sc=15){
  types_int <- c('BP', 'MF', 'CC')
  onts <- c('Biological Process', 'Molecular Function', 'Cellular Component')
  side <- c('Right', 'Left')
  plot <- enrichGO(gene = genelist, 
           OrgDb = org.EcK12.eg.db,
           keyType = "ENTREZID", 
           ont = types_int[ont_type],
           pvalueCutoff = pcut_gp,
           pAdjustMethod = "BH",
           qvalueCutoff = qcut_gp,
           readable = TRUE) %>% goplot(showCategory = sc) +
    ggtitle(paste0(side[side_num], ' genes, ',onts[ont_type] ,' onthology, nd=', as.character(nd)))+
    theme(plot.title = element_text(hjust = 0.5))
  return(plot)
}

#nd - neighbor distance
get_plots <- function(nd, pcut=1, qcut=1){
  neighbors <- get_neighbor(nd)
  go_right <- neighbors %>% dplyr::filter(side=='right') %>% pull(ENTREZID) %>% na.omit() %>% as.vector()
  go_left <- neighbors %>% dplyr::filter(side=='left') %>% pull(ENTREZID) %>% na.omit()%>% as.vector()
  ego_right <- enrichGO(gene = go_right, 
                        OrgDb = org.EcK12.eg.db,
                        keyType = "ENTREZID", 
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
  right_go_bp <- get_goplot(go_right, 1, 1, nd, pcut_gp=pcut)
  right_go_mf <- get_goplot(go_right, 2, 1, nd, pcut_gp=pcut)
  right_go_cc <- get_goplot(go_right, 3, 1, nd, pcut_gp=pcut)
  
  ego_left <- enrichGO(gene = go_left, 
                        OrgDb = org.EcK12.eg.db,
                        keyType = "ENTREZID", 
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
  left_go_bp <- get_goplot(go_left, 1, 2, nd, pcut_gp=pcut)
  left_go_mf <- get_goplot(go_left, 2, 2, nd, pcut_gp=pcut)
  left_go_cc <- get_goplot(go_left, 3, 2, nd, pcut_gp=pcut)
  
  right_kegg <- neighbors %>% dplyr::filter(side=='right')%>% pull(kegg)%>% na.omit() %>% as.vector()
  left_kegg <- neighbors %>% dplyr::filter(side=='left')%>% pull(kegg)%>% na.omit() %>% as.vector()
  ekegg_right <- enrichKEGG(gene = right_kegg,
                            organism     = 'eco',
                            pvalueCutoff = 1)
  ekegg_left <- enrichKEGG(gene = left_kegg,
                            organism     = 'eco',
                            pvalueCutoff = 1)
  if (nrow(ekegg_right@result[ekegg_right@result$p.adjust<=0.05,]) == 0) {
    dot_kegg_right <- ggplot() + 
      ggtitle(paste0("No KEGG enrichment found for right genes, nd=", nd))
  } else {
    dot_kegg_right <- dotplot(ekegg_right, title = paste0('Right genes, KEGG, nd=', nd))
  }
  
  if (nrow(ekegg_left@result[ekegg_left@result$p.adjust<=0.05,]) == 0) {
    dot_kegg_left <- ggplot() + 
      ggtitle(paste0("No KEGG enrichment found for left genes, nd=", nd))
  } else {
    dot_kegg_left <- dotplot(ekegg_left, title = paste0('Left genes, KEGG, nd=', nd))
  }
  
  return(list(right_dots, right_go_mf, right_go_cc, right_go_bp, left_dots, left_go_mf, left_go_cc, left_go_bp,dot_kegg_right, dot_kegg_left))
}

pdf("GO_enrich_plots_result.pdf", width = 15, height = 12)
get_plots(nd=1, pcut=0.5)
get_plots(nd=2, pcut=0.3)
get_plots(nd=3, pcut=0.3)
get_plots(nd=5, pcut=0.2)
dev.off()



#View
neighbors <- get_neighbor(nd=3)
go_right <- neighbors %>% dplyr::filter(side=='right') %>% pull(ENTREZID) %>% na.omit()%>% as.vector()
go_left <- neighbors %>% dplyr::filter(side=='left') %>% pull(ENTREZID) %>% na.omit() %>% as.vector()
ego_right <- enrichGO(gene = go_right, 
                      OrgDb = org.EcK12.eg.db,
                      keyType = "ENTREZID", 
                      ont = "ALL",
                      pvalueCutoff = 1,
                      pAdjustMethod = "BH",
                      qvalueCutoff = 1,
                      readable = TRUE)
ego_left <- enrichGO(gene = go_left, 
                      OrgDb = org.EcK12.eg.db,
                      keyType = "ENTREZID", 
                      ont = "ALL",
                      pvalueCutoff = 1,
                      pAdjustMethod = "BH",
                      qvalueCutoff = 1,
                      readable = TRUE)

ego_results <- rbind(ego_right@result %>% mutate(side='right'), ego_left@result %>% mutate(side='left')) %>% 
  filter(p.adjust<=0.05) %>% arrange(desc(Count))

#top genes
ego_mf <- ego_results %>% filter(ONTOLOGY=='MF') %>% filter(Count > quantile(Count, probs = seq(0, 1, 0.1))[10])
ego_bp <- ego_results %>% filter(ONTOLOGY=='BP') %>% filter(Count > quantile(Count, probs = seq(0, 1, 0.1))[10])
ego_cc <- ego_results %>% filter(ONTOLOGY=='CC') %>% filter(Count > quantile(Count, probs = seq(0, 1, 0.1))[10])

#propionate
propionate1 <- ego_right@result %>% dplyr::filter(str_detect(Description, "propionate")) %>% filter(p.adjust<=0.05) %>% mutate(side='right')
propionate2 <- ego_left@result %>% dplyr::filter(str_detect(Description, "propionate")) %>% filter(p.adjust<=0.05) %>% mutate(side='left')
rbind(propionate1, propionate2) %>% arrange(p.adjust) #%>% write_tsv('propionate.tsv')

#kegg res
right_kegg <- neighbors %>% dplyr::filter(side=='right')%>% pull(kegg)%>% na.omit() %>% as.vector()
left_kegg <- neighbors %>% dplyr::filter(side=='left')%>% pull(kegg)%>% na.omit() %>% as.vector()
ekegg_right <- enrichKEGG(gene = right_kegg,
                          organism = 'eco',
                          pvalueCutoff = 1)
ekegg_left <- enrichKEGG(gene = left_kegg,
                         organism     = 'eco',
                         pvalueCutoff = 1)
kegg_res <- rbind(ekegg_right@result %>% mutate(side='right'), ekegg_left@result %>% mutate(side='left')) %>% 
  filter(p.adjust<=0.05) %>% arrange(desc(Count))

#web pics
for (ids  in (kegg_res %>% filter(side=='left') %>% pull(ID))){
  browseKEGG(ekegg_left, ids)
}
for (ids  in (kegg_res %>% filter(side=='right') %>% pull(ID))){
  browseKEGG(ekegg_right, ids)
}

#strsplit
neighbors %>% pull(product) %>% strsplit(., ' ') %>% unlist() %>% strsplit(., '-') %>% unlist() %>% table() %>% sort()

words <- neighbors %>% pull(product) %>% strsplit(., ' ') %>% unlist() %>% strsplit(., '-') %>% 
  unlist() %>% table() %>% sort(decreasing = T) %>% as.data.frame() %>% filter(Freq>1) %>% 
  rename('.' = 'wors') %>% mutate(str_len=str_length(wors)) %>% filter(str_len>=3)
