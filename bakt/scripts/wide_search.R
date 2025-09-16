library(tidyverse)
library(clusterProfiler)
library(DOSE)
library(enrichplot)
library(org.EcK12.eg.db)

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

get_neighbor <- function(nd, x){
  idxes <- which(x$gene == 'naRNA4')
  right_neighbors <- x %>% dplyr::slice(outer((1:nd),idxes, '+') %>% as.vector()) %>% dplyr::filter(gene!='naRNA4') %>% mutate(side='right')
  left_neighbors <- x %>% dplyr::slice(outer((-nd:-1),idxes, '+') %>% as.vector()) %>% dplyr::filter(gene!='naRNA4')%>% mutate(side='left')
  neighbors <- rbind(right_neighbors, left_neighbors) %>% dplyr::arrange(start) %>% distinct(start, .keep_all = TRUE)
  return(neighbors)
}

tab_analysis <- function(x, n){
  neig <- get_neighbor(n, x)
  go_left <- neig %>% filter(side=='left') %>% pull(ENTREZID)%>% na.omit() %>% as.character()
  go_right <- neig %>% filter(side=='right') %>% pull(ENTREZID)%>% na.omit() %>% as.character()
  ego_right <- enrichGO(gene = go_right, 
                        OrgDb = org.EcK12.eg.db,
                        keyType = "ENTREZID", 
                        ont = "ALL",
                        pvalueCutoff = 0.05,
                        pAdjustMethod = "BH",
                        qvalueCutoff = 1,
                        readable = TRUE)
  ego_left <- enrichGO(gene = go_right, 
                       OrgDb = org.EcK12.eg.db,
                       keyType = "ENTREZID", 
                       ont = "ALL",
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "BH",
                       qvalueCutoff = 1,
                       readable = TRUE)
  right_kegg <- neig %>% dplyr::filter(side=='right')%>% pull(kegg)%>% na.omit() %>% as.character()
  left_kegg <- neig %>% dplyr::filter(side=='left')%>% pull(kegg)%>% na.omit() %>% as.character()
  ekegg_right <- enrichKEGG(gene = right_kegg,
                            organism     = 'eco',
                            pvalueCutoff = 0.05)
  ekegg_left <- enrichKEGG(gene = left_kegg,
                           organism     = 'eco',
                           pvalueCutoff = 0.05)
  return(c(ego_right, ego_left, ekegg_right, ekegg_left))
}

get_go_plots <- function(row_n, title_char, list_of_res){
  tablet <- list_of_res[[row_n,]] %>% 
    map(~ .x@result) %>%           # извлекаем @result из каждого элемента
    imap(~ mutate(.x, file_name = .y)) %>%  # добавляем имя файла
    bind_rows() %>% filter(ONTOLOGY == 'BP')%>% mutate(file_name =   sub('.*/', '', file_name)) %>% 
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

tablo <- readRDS("all_narna.Rds")
imp <- c('ECZV_18330','ECZV_18305','ECZV_18395','ECZV_20640','ECZV_19270')
targets <- tablo %>% filter(seq_name %in% paste(imp, 'Nucleoid-associated noncoding RNA 4 (CssrE)')) %>% pull(new_name)

bakts <- tablo %>% filter(new_name %in% targets) %>% pull(seq_name) %>% gsub(' Nucleoid-associated noncoding RNA 4 \\(CssrE\\)', '', .) %>% gsub('_\\d+$', '', .) %>%  unique()
#tablo %>% filter(new_name %in% targets) %>% pull(seq_name) %>%  gsub(' Nucleoid-associated noncoding RNA 4 \\(CssrE\\)', '', .) %>% gsub('_\\d+$', '', .) %>% table()

contain_dir <- dir('bakta_isolates')[dir('bakta_isolates') %in% paste0('bakta_annotation_',bakts)]
contained_tsv <- list.files(file.path('bakta_isolates', contain_dir), full.names = TRUE) %>% grep('\\d+\\.tsv$', ., value = T)
gen_tables <- lapply(contained_tsv, read_anat_tab)

tab_res_nd1 <- lapply(gen_tables, function(x){tab_analysis(x, 1)})
sapply(tab_res_nd1, function(x){get_go_plots(1, 'test', x)})

