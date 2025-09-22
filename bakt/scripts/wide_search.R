library(tidyverse)
library(clusterProfiler)
library(DOSE)
library(enrichplot)
library(org.EcK12.eg.db)

#Задача функций----------------------
read_anat_tab <- function(path){
  lines <- readLines(path)
  hash_lines <- grep("^#", lines)
  skip_lines <- if (length(hash_lines) > 0) max(hash_lines) - 1 else 0
  eczv <- read_tsv(path, skip = skip_lines) %>% janitor::clean_names()
  file_name <- sub('.*/', '', path)
  entr <- bitr(
    eczv$gene,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.EcK12.eg.db) %>% dplyr::rename('gene' = 'SYMBOL')
  eczv <- left_join(eczv, entr, by='gene')
  keggs <- bitr_kegg(eczv$ENTREZID, fromType='ncbi-geneid', toType='kegg', organism='eco', drop = TRUE) %>% dplyr::rename('ENTREZID' = 'ncbi-geneid' )
  eczv <- left_join(eczv, keggs, by='ENTREZID')
  return(eczv)
}

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
        li <- li+1
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
        ri <- ri+1
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

get_enrichment <- function(gene_set, db_type='GO', pcut=0.05, qcut=0.05){
  db_types <- c('GO', 'KEGG')
  if (!(db_type %in% db_types)){
    print('Not available enrichment type')
    return(NA)
  }
  if (db_type=='GO'){
    enrich <- enrichGO(gene = gene_set, 
             OrgDb = org.EcK12.eg.db,
             keyType = "ENTREZID", 
             ont = "ALL",
             pvalueCutoff = pcut,
             pAdjustMethod = "BH",
             qvalueCutoff = qcut,
             readable = TRUE)
  } else if(db_type=='KEGG'){
    enrich <- enrichKEGG(gene = gene_set,
               organism     = 'eco',
               pvalueCutoff = pcut)
  }
  return(enrich)
}

enrich_to_db_side_tab <- function(x, side_str, ontology_str){
  if (is.null(x) || !isS4(x) || is.null(slot(x, "result"))) {
    return(NULL)
  } else{x@result %>% mutate(side=side_str, ontology=ontology_str) %>% return()}
}

tab_analysis <- function(x, n, loc_tag='naRNA4', pcut=0.05, qcut=0.05, alt_tab=T){
  if (alt_tab) {
    neig <- alt_neighb(n, x, tag=loc_tag)
  } else{
    neig <- get_neighbor(n, x, tag=loc_tag)
  }
  go_left <- neig %>% filter(side=='left') %>% pull(ENTREZID) %>% na.omit() %>% as.character()
  go_right <- neig %>% filter(side=='right') %>% pull(ENTREZID) %>% na.omit() %>% as.character()
  right_kegg <- neig %>% dplyr::filter(side=='right')%>% pull(kegg)%>% na.omit() %>% as.character()
  left_kegg <- neig %>% dplyr::filter(side=='left')%>% pull(kegg)%>% na.omit() %>% as.character()
  
  print(head(go_right))
  ego_right <- get_enrichment(go_right, db_type='GO', pcut=0.05, qcut=0.05) %>% enrich_to_db_side_tab(side_str = 'right', ontology_str = 'GO')
  print(head(go_left))
  ego_left <- get_enrichment(go_left, db_type='GO', pcut=0.05, qcut=0.05) %>% enrich_to_db_side_tab(side_str = 'left', ontology_str = 'GO')
  print(head(right_kegg))
  ekegg_right <- get_enrichment(right_kegg, db_type='KEGG', pcut=0.05, qcut=0.05) %>% enrich_to_db_side_tab(side_str = 'right', ontology_str = 'KEGG')
  print(head(left_kegg))
  ekegg_left <- get_enrichment(left_kegg, db_type='KEGG', pcut=0.05, qcut=0.05) %>% enrich_to_db_side_tab(side_str = 'left', ontology_str = 'KEGG')

  en_res_list <- list(ego_right, ego_left, ekegg_right, ekegg_left) %>% 
    keep(~ !is.null(.x)) %>% keep(~ is.data.frame(.x) && (nrow(.x) > 0))
  
  if (length(en_res_list) == 0) {
    message("--> All enrichment results are empty or NULL. Returning empty data.frame.")
    return(data.frame())
  }
  
  en_res <- bind_rows(en_res_list)
  return(en_res)
}

#Предобработка-----------------

isols <- read_tsv("AIEC_isol.tsv") #Статусы патогенности изолятов
tablo <- readRDS("all_narna.Rds") #Обобщённая аннотация
imp <- c('ECZV_18330','ECZV_18305','ECZV_18395','ECZV_20640','ECZV_19270') # 5 самых значимых генов
targets <- tablo %>% filter(seq_name %in% paste(imp, 'Nucleoid-associated noncoding RNA 4 (CssrE)')) %>% pull(new_name) #получение обобщённых названий генов среди всех изолятов

target_seqs <- tablo %>% filter(new_name %in% targets) %>% pull(seq_name) %>% gsub(' Nucleoid-associated noncoding RNA 4 \\(CssrE\\)', '', .) #Извлечение id генов внутри изолятов
bakts <- tablo %>% filter(new_name %in% targets) %>% pull(seq_name) %>% gsub(' Nucleoid-associated noncoding RNA 4 \\(CssrE\\)', '', .) %>% gsub('_\\d+$', '', .) %>%  unique() #Получение списка изолятова
#tablo %>% filter(new_name %in% targets) %>% pull(seq_name) %>%  gsub(' Nucleoid-associated noncoding RNA 4 \\(CssrE\\)', '', .) %>% gsub('_\\d+$', '', .) %>% table()

contain_dir <- dir('bakta_isolates')[dir('bakta_isolates') %in% paste0('bakta_annotation_',bakts)]
contained_tsv <- list.files(file.path('bakta_isolates', contain_dir), full.names = TRUE) %>% grep('\\d+\\.tsv$', ., value = T)
gen_tables <- lapply(contained_tsv, read_anat_tab)


#Обработка-----------

fultab_res_nd1 <- lapply(gen_tables, function(x){tab_analysis(x, n=1, loc_tag=target_seqs, pcut=0.1, qcut=0.1)})
names(fultab_res_nd1) <- basename(contained_tsv) %>% sub('.tsv', '', ., fixed = T)
sapply(fultab_res_nd1, function(x){if (nrow(x) > 1) {
  return(x %>% pull(Count))
} else {
  return(NA)}
}) %>% unlist() %>% hist()
sapply(fultab_res_nd1, function(x){if (nrow(x) > 1) {
  return(x %>% pull(geneID) %>% strsplit('/') %>% unlist())
} else {
  return(NA)}
}) %>% unlist() %>% table() %>% sort()

fultab_res_nd3 <- lapply(gen_tables, function(x){tab_analysis(x, n=3, loc_tag=target_seqs, pcut=0.1, qcut=0.1)})
names(fultab_res_nd3) <- basename(contained_tsv) %>% sub('.tsv', '', ., fixed = T)
sapply(fultab_res_nd3, function(x){if (nrow(x) > 1) {
  return(x %>% pull(Count))
} else {
  return(NA)}
}) %>% unlist() %>% hist()
sapply(fultab_res_nd3, function(x){if (nrow(x) > 1) {
  return(x %>% pull(geneID) %>% strsplit('/') %>% unlist())
} else {
  return(NA)}
}) %>% unlist() %>% table() %>% sort()
sapply(fultab_res_nd3, function(x){
  x %>% filter(Count>2 & p.adjust<=0.05)
}) %>% bind_rows() %>% pull(Description) %>% table()

plus_tabs_nd3 <- lapply(gen_tables, function(x){tab_analysis((x %>% filter(strand=='+')), n=3, loc_tag=target_seqs, pcut=0.1, qcut=0.1)})
names(plus_tabs_nd3) <- basename(contained_tsv) %>% sub('.tsv', '', ., fixed = T)
sapply(plus_tabs_nd3, function(x){
  if(nrow(x)>0){x  %>% filter(Count>1 & p.adjust<=0.01)}}) %>% 
  bind_rows() %>% pull(Description) %>% table()%>% sort(decreasing = T)
sapply(plus_tabs_nd3, function(x){
  if(nrow(x)>0){x  %>% filter(Count>1 & p.adjust<=0.01)}}) %>% 
  bind_rows() %>% pull(side) %>% table()

minus_tabs_nd3 <- lapply(gen_tables, function(x){tab_analysis((x %>% filter(strand=='-')), n=3, loc_tag=target_seqs, pcut=0.1, qcut=0.1)})
names(minus_tabs_nd3) <- basename(contained_tsv) %>% sub('.tsv', '', ., fixed = T)
sapply(minus_tabs_nd3, function(x){
  if(nrow(x)>0){x  %>% filter(Count>1 & p.adjust<=0.01)}}) %>% 
  bind_rows() %>% pull(Description) %>% table() %>% sort(decreasing = T)
sapply(minus_tabs_nd3, function(x){
  if(nrow(x)>0){x  %>% filter(Count>1 & p.adjust<=0.01)}}) %>% 
  bind_rows() %>% pull(side) %>% table()

#sapply(plus_tabs_nd3, nrow) %>% as.data.frame() %>% rownames_to_column(var='File') %>% rename('.'= 'n_row') %>% left_join(isols, by='File') %>% dplyr::select(1, 2, 5) %>% dplyr::slice(-1) 
#Имеет смысл сделать это раньше

minus_tabs_nd3 <- lapply(gen_tables, function(x){tab_analysis((x %>% filter(strand=='-')), n=1, loc_tag=target_seqs, pcut=0.1, qcut=0.1)})
names(minus_tabs_nd1) <- basename(contained_tsv) %>% sub('.tsv', '', ., fixed = T)


#garb------
get_go_plots <- function(row_n, title_char, list_of_res){
  tablet <- list_of_res[row_n] %>% 
    map(~ .x@result) %>%
    # Create file_name and process it in the same pipeline step
    imap(~ mutate(.x, file_name = sub('.*/', '', .y))) %>%  # Process file name immediately
    bind_rows() %>%
    # Now file_name exists and can be further processed
    mutate(ratio = as.numeric(sub("/.*", "", GeneRatio)) / as.numeric(sub(".*/", "", GeneRatio)))
  plot <- ggplot()+
    geom_point(data=tablet, aes(x=ratio, y=Description, color=p.adjust, size=Count))+
    #geom_text(data=tablet, aes(x=ratio, y=Description, label=file_name), 
    #size=2, hjust=0, vjust=0.1, nudge_x=0.02)+
    scale_color_gradientn(colors = c("darkblue", "red"))+
    scale_size_continuous(range = c(5, 8))+
    labs(title=title_char)+
    xlab('Gene ratio')
  return(plot)
}
get_kegg_plots <- function(row_n, title_char, list_of_res){
  tablet <- list_of_res[row_n,] %>% 
    map(~ .x@result) %>%
    imap(~ mutate(.x, file_name = .y)) %>%
    bind_rows() %>% filter(p.adjust<=0.05) %>% mutate(file_name =   sub('.*/', '', file_name)) %>% 
    mutate(ratio = as.numeric(sub("/.*", "", GeneRatio)) / as.numeric(sub(".*/", "", GeneRatio)))
  plot <- ggplot()+
    geom_point(data=tablet, aes(x=ratio, y=Description, color=p.adjust, size=Count))+
    #geom_text(data=tablet, aes(x=ratio, y=Description, label=file_name), 
    #size=2, hjust=0, vjust=0.1, nudge_x=0.02)+
    scale_color_gradientn(colors = c("darkblue", "red"))+
    scale_size_continuous(range = c(5, 8))+
    labs(title=title_char)+
    xlab('Gene ratio')
  return(plot)
}