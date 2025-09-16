library(tidyverse)
library(reshape2)
library(ComplexHeatmap)
library(dendextend)

all_names <- readRDS('all_narna.Rds')
text <- readLines("all_naRNA_uniq_dist.mat")
text <- text[-c(grep("^>", text) %>% max(), length(text))]

names <- sub("^>", "", grep("^>", text, value = TRUE)) #seq names
#names <- data.frame('seq_name'= names) %>% left_join(all_names, by='seq_name') %>% dplyr::pull(new_name)
names <- sub('\\sNucleoid-associated noncoding RNA 4 \\(CssrE\\)$', '', names)

data_lines <- text[!grepl("^>", text) & text != ""] #dist values
split_lines <- lapply(data_lines, function(line) as.numeric(strsplit(trimws(line), "\\s+")[[1]]))
dist_values <- unlist(split_lines)
n <- length(split_lines) + 1

dist_matrix <- matrix(0, nrow = n, ncol = n, byrow = TRUE)
idx <- which(lower.tri(dist_matrix), arr.ind = TRUE)
idx_ordered <- idx[order(idx[,1], idx[,2]), ]

for (i in seq_len(nrow(idx_ordered))) {
  dist_matrix[idx_ordered[i, "row"], idx_ordered[i, "col"]] <- dist_values[i]
}

dist_matrix <- dist_matrix + t(dist_matrix)
rownames(dist_matrix) <- colnames(dist_matrix) <- names


#heatmap
Heatmap(dist_matrix, row_names_gp = gpar(fontsize = 5),column_names_gp = gpar(fontsize = 5), heatmap_legend_param = list(title = 'base pairs \n distance'))
heatmap(dist_matrix, symm = TRUE, main = "Distance Matrix Heatmap", cexRow = 0.5, cexCol = 0.5)

#Дерево
get_color <- function(name, col1='red', col2='blue') {
  if (grepl("-\\d{5}$", name)) {
    return(col1)
  } else if (grepl("_\\d{5}$", name)) {
    return(col2)
  } else {
    return("black")
  }
}

hc <- hclust(as.dist(dist_matrix))
hc$height <- hc$height / max(hc$height)
dend <- as.dendrogram(hc)
labels_colors(dend) <- sapply(labels(dend), get_color)
plot(dend, main = "Hierarchical Clustering with Colored Labels")
