library(tidyverse)
library(clusterProfiler)
library(DOSE)
library(enrichplot)
library(org.EcK12.eg.db)
library(readxl)

eczv <- read_tsv('bakt_files/ECZV.tsv', skip=2) %>% janitor::clean_names()
rna_seq <- read_excel("rna_dif_exp.xlsx",sheet=1) %>% rename('gene' = 'locus_tag') %>% mutate(across(6:11, as.double))
eczv <- left_join(eczv, rna_seq, by='locus_tag')

entr <- bitr(
  eczv$gene,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.EcK12.eg.db) %>% dplyr::rename('gene' = 'SYMBOL')
eczv <- left_join(eczv, entr, by='gene')
keggs <- bitr_kegg(eczv$ENTREZID, fromType='ncbi-geneid', toType='kegg', organism='eco', drop = TRUE) %>% dplyr::rename('ENTREZID' = 'ncbi-geneid' )
eczv <- left_join(eczv, keggs, by='ENTREZID')

get_neighbor <- function(nd, x, tag='naRNA4'){
  if (any(tag == 'naRNA4')) {
    idxes <- which(x$gene == 'naRNA4')
  } else{
    idxes <- which(x$locus_tag %in% tag)
  }
  plus_strand <- x %>% filter(strand=='+')
  minus_strand <- x %>% filter(strand=='-')
  plus_right_neighbors <- plus_strand %>% dplyr::slice(outer((1:nd),idxes, '+') %>% as.vector()) %>% dplyr::filter(gene!='naRNA4') %>% mutate(side='right')
  plus_left_neighbors <- plus_strand %>% dplyr::slice(outer((-nd:-1),idxes, '+') %>% as.vector()) %>% dplyr::filter(gene!='naRNA4')%>% mutate(side='left')
  minus_right_neighbors <- minus_strand %>% dplyr::slice(outer((1:nd),idxes, '+') %>% as.vector()) %>% dplyr::filter(gene!='naRNA4') %>% mutate(side='right')
  minus_left_neighbors <- minus_strand %>% dplyr::slice(outer((-nd:-1),idxes, '+') %>% as.vector()) %>% dplyr::filter(gene!='naRNA4')%>% mutate(side='left')
  neighbors <- rbind(plus_right_neighbors, plus_left_neighbors, minus_right_neighbors, minus_left_neighbors) %>% dplyr::arrange(start) %>% distinct(start, .keep_all = TRUE)
  return(neighbors)
}

alt_neighb <- function(nd, x, tag='naRNA4') {
  tab <- x %>% mutate(rnum=row_number()) %>% filter(!is.na(gene))
  if (any(tag == 'naRNA4')) {
    idxes <- which(tab$gene == 'naRNA4')
  } else{
    idxes <- which(tab$locus_tag %in% tag)
  }
  left_neighb <- list()
  right_neighb <- list()
  for (id in idxes){
    left <- integer()
    right <- integer()
    li <- 1
    ri <- 1
    while (length(left)!=nd){
      #if (id-li > nrow(x)) {
      #  break
      #}
      if (tab$gene[id-li]!='naRNA4'){
        left <- base::append(left, id-li)
      } else{
        li <- li+1
      }
    }
    while (length(right)!=nd){
      #if (id+ri > nrow(x)) {
      #  break  # или stop("Недостаточно строк справа")
      #}
      if (tab$gene[id+ri]!='naRNA4'){
        right <- base::append(right, id+ri)
      } else{
        ri <- ri+1
      }
    }
    left_neighb <- base::append(left_neighb, list(left))
    right_neighb <- base::append(right_neighb, list(right))
  }
  tab_right <- tab %>% dplyr::slice(unique(unlist(right_neighb))) %>% mutate(side='right')
  tab_left <- tab %>% dplyr::slice(unique(unlist(left_neighb)))%>% mutate(side='left')
  return(bind_rows(tab_right, tab_left))
}

get_plots <- function(nd, x, alt=T, loc_tag='naRNA4', title=''){
  neib <- if(alt){alt_neighb(1, x, tag=loc_tag)}else{get_neighbor(nd, x, tag=loc_tag)}
  fc_vec <-  neib %>%  group_by(ENTREZID) %>% 
    summarise(log2FoldChange = median(log2FoldChange)) %>% 
    dplyr::select(ENTREZID, log2FoldChange) %>% 
    na.omit() %>% deframe()  %>% sort(decreasing = T)
  gsea_go <-  gseGO(gene=fc_vec,
          OrgDb= org.EcK12.eg.db,
          keyType= 'ENTREZID',
          ont= "BP")
  gsea_kegg <-  neib %>%  group_by(ENTREZID) %>% 
    summarise(log2FoldChange = median(log2FoldChange)) %>% 
    dplyr::select(ENTREZID, log2FoldChange) %>% 
    na.omit() %>% deframe()  %>% sort(decreasing = T) %>% 
    gseKEGG(geneList=.,
            organism='eco',
            keyType='ncbi-geneid',
            pvalueCutoff = 0.05,
            verbose= FALSE)
    gop <-  setReadable(pairwise_termsim(gsea_go), 'org.EcK12.eg.db', 'ENTREZID')
  dog <- dotplot(gsea_go)+
    ggtitle(paste0('GO ontology, nd=', nd, title))+
    theme(plot.title = element_text(hjust = 0.5))
  gog <- goplot(gsea_go)+
    ggtitle(paste0('GO ontology, nd=', nd, title))
  goh <- heatplot(gop, showCategory=10, foldChange = fc_vec)+
    ggtitle(paste0('GO ontology, nd=', nd, title))+
    theme(plot.title = element_text(hjust = 0.5))
  goc <- cnetplot(gop, showCategory=5, foldChange = fc_vec,
           cex_label_gene = 0.5)+
    ggtitle(paste0('GO ontology, nd=', nd, title))+
    theme(plot.title = element_text(hjust = 0.5))
  got <- treeplot(gop)+
    ggtitle(paste0('GO ontology, nd=', nd, title))+
    theme(plot.title = element_text(hjust = 0.5))
  
    kep <- pairwise_termsim(gsea_kegg)
  ked <- dotplot(gsea_kegg)+
    ggtitle(paste0('KEGG ontology, nd=', nd, title))+
    theme(plot.title = element_text(hjust = 0.5))
  
  sapply(list(dog, gog, goh, goc, got, ked), print)
}

#get_plots(1, eczv, title=' cluster`s edges')
get_plots(3, eczv, title=', ALL naRNA4, cluster`s edges')
get_plots(3, eczv %>% filter(strand=='-'), title=', ALL naRNA4, cluster`s edges, minus strand')
