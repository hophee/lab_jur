library(tidyverse)
library(clusterProfiler)
library(org.EcK12.eg.db)

#hand made
is_canon <- c(0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1)

neighbor_of_canon <- function(nd, tab=eczv, old=F){
  idx_canon <- which((tab$gene == 'naRNA4')&(tab$is_canon==1))
  idx_non_canon <- which((tab$gene == 'naRNA4')&(tab$is_canon==0))
  can <- tab %>% dplyr::slice(c(outer((1:nd),idx_canon, '+'), outer((-nd:-1),idx_canon, '+')) %>% as.vector()) %>% mutate(canon_neighb=1)
  no_can <- tab %>% dplyr::slice(c(outer((1:nd),idx_non_canon, '+'), outer((-nd:-1),idx_non_canon, '+')) %>% as.vector()) %>% mutate(canon_neighb=0)
  if(old){
    neighbors <- rbind(can, no_can) %>% dplyr::arrange(start)
  }else{
    neighbors <- rbind(can, no_can) %>% dplyr::arrange(start) %>% distinct(start, .keep_all = TRUE) %>% filter(gene!='naRNA4')
  }
  return(neighbors)
}

eczv <- read_tsv('bakt_files/ECZV.tsv', skip=2) %>% janitor::clean_names()
locus_tag <- eczv %>% filter(gene=='naRNA4') %>% pull(locus_tag) %>% sort()
canon_tab <- data.frame(locus_tag, is_canon) #%>% dplyr::slice(-58) #в 58 строке спорный образец
old_all_neighb <- left_join(eczv, canon_tab, by='locus_tag') %>% neighbor_of_canon(1, tab=., old=T) %>% dplyr::select(start, type, canon_neighb, gene, locus_tag)

non_can <- old_all_neighb %>% filter(type %in% c('cds', 'ncRNA')) %>% filter(canon_neighb==0) %>% pull(type) %>% as.factor()
can <- old_all_neighb %>% filter(type %in% c('cds', 'ncRNA')) %>% filter(canon_neighb==1) %>% pull(type) %>% as.factor()
contingency_table <- table(group = c(rep("non_can", length(non_can)), rep("can", length(can))),
                           value = c(non_can, can))
chi2_test <- chisq.test(contingency_table)
#Значение p-value
print(chi2_test$p.value)
#Наблюдаемые значения
print(chi2_test$observed)

#---------
all_neighb <- left_join(eczv, canon_tab, by='locus_tag') %>% neighbor_of_canon(1, tab=., old=F) %>% dplyr::select(start, type, canon_neighb, gene, locus_tag)
all_neighb_entr <- all_neighb %>%
  left_join(
    bitr(all_neighb$gene, 
         fromType = "SYMBOL", 
         toType = "ENTREZID", 
         OrgDb = org.EcK12.eg.db),
    by = c("gene" = "SYMBOL")
  ) %>%
  na.omit()
canon_genes <- all_neighb_entr %>% 
  filter(canon_neighb == 1) %>% 
  pull(ENTREZID) %>% 
  unique()
non_canon_genes <- all_neighb_entr %>% 
  filter(canon_neighb == 0) %>% 
  pull(ENTREZID) %>% 
  unique()
gene_clusters <- list(
  "Canonical" = canon_genes,
  "Non-canonical" = non_canon_genes
)


all_go_terms <- compareCluster(
  geneCluster = gene_clusters,
  fun = "enrichGO",
  OrgDb = org.EcK12.eg.db,
  keyType = "ENTREZID",
  ont = 'ALL',
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH")
dotplot(all_go_terms, showCategory = 15, font.size = 10) +
  ggtitle("GO Terms: Canonical vs Non-canonical neighbors")

##Splited same
# clusters <- sapply(c('BP', 'MF'), function(x){compareCluster(
#   geneCluster = gene_clusters,
#   fun = "enrichGO",
#   OrgDb = org.EcK12.eg.db,
#   keyType = "ENTREZID",
#   ont = x,
#   pvalueCutoff = 0.05,
#   pAdjustMethod = "BH"
# )})
# dotplot(clusters[[1]], showCategory = 15, font.size = 10) +
#   ggtitle("GO BP: Canonical vs Non-canonical neighbors")
# dotplot(clusters[[2]], showCategory = 15, font.size = 10) +
#   ggtitle("GO MF: Canonical vs Non-canonical neighbors")


#KEGG обогащение
kegg_mapping <- bitr_kegg(
  gene = unlist(gene_clusters),
  fromType = 'ncbi-geneid',
  toType = 'kegg',
  organism = 'eco',
  drop = TRUE
)
# Разделяем KEGG ID по группам
kegg_clusters <- lapply(gene_clusters, function(entrez_ids) {
  # Находим KEGG ID для генов из этой группы
  kegg_ids <- kegg_mapping %>%
    filter(`ncbi-geneid` %in% entrez_ids) %>%
    pull(kegg) %>%
    unique()  # Убираем дубликаты
  return(kegg_ids)
})
kegg_comp <- compareCluster(
  geneCluster = kegg_clusters,
  fun = "enrichKEGG",
  organism = "eco",
  pvalueCutoff = 0.05
)
dotplot(kegg_comp, showCategory = 15, font.size = 10) +
  ggtitle("KEGG: Canonical vs Non-canonical neighbors")

#-------------
#Same, but split strands

plots_of_enrich <- function(tab, add_text=''){
  neighb <- left_join(tab, canon_tab, by='locus_tag') %>% neighbor_of_canon(1, tab=.) %>% dplyr::select(start, type, canon_neighb, gene, locus_tag)
  neighb_entr <- neighb %>%
    left_join(
      bitr(neighb$gene, 
           fromType = "SYMBOL", 
           toType = "ENTREZID", 
           OrgDb = org.EcK12.eg.db),
      by = c("gene" = "SYMBOL")
    ) %>%
    na.omit()
  canon_genes <- neighb_entr %>% 
    filter(canon_neighb == 1) %>% 
    pull(ENTREZID) %>% 
    unique()
  non_canon_genes <- neighb_entr %>% 
    filter(canon_neighb == 0) %>% 
    pull(ENTREZID) %>% 
    unique()
  gene_clusters <- list(
    "Canonical" = canon_genes,
    "Non-canonical" = non_canon_genes
  )
  
  print('Compairing clusters, GO terms')
  all_go_terms <- compareCluster(
    geneCluster = gene_clusters,
    fun = "enrichGO",
    OrgDb = org.EcK12.eg.db,
    keyType = "ENTREZID",
    ont = 'ALL',
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH")
  if (!is.null(all_go_terms)){dot1 <- dotplot(all_go_terms, showCategory = 15, font.size = 10) +
    ggtitle(paste0("GO Terms: Canonical vs Non-canonical neighbors", add_text))} else{
      dot1 <- 'No GO enrichment'
    }
  
  
  kegg_mapping <- bitr_kegg(
    gene = unlist(gene_clusters),
    fromType = 'ncbi-geneid',
    toType = 'kegg',
    organism = 'eco',
    drop = TRUE
  )
  # Разделяем KEGG ID по группам
  kegg_clusters <- lapply(gene_clusters, function(entrez_ids) {
    # Находим KEGG ID для генов из этой группы
    kegg_ids <- kegg_mapping %>%
      filter(`ncbi-geneid` %in% entrez_ids) %>%
      pull(kegg) %>%
      unique()  # Убираем дубликаты
    return(kegg_ids)
  })
  print('Compairing clusters, KEGG')
  kegg_comp <- compareCluster(
    geneCluster = kegg_clusters,
    fun = "enrichKEGG",
    organism = "eco",
    pvalueCutoff = 0.05
  )
  if (!is.null(kegg_comp)){dot2 <- dotplot(kegg_comp, showCategory = 15, font.size = 10) +
    ggtitle(paste0("KEGG Terms: Canonical vs Non-canonical neighbors", add_text))}else{
      dot2 <- 'No KEGG enrichment'
    }
  
  print(dot1)
  print(dot2)
}

eczv_right <- eczv %>% filter(strand=='+')
eczv_left <- eczv %>% filter(strand=='-')

plots_of_enrich(eczv_right, add_text = '; plus strand')
plots_of_enrich(eczv_left, add_text = '; minus strand')

