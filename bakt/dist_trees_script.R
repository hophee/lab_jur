library(tidyverse)
library(dendextend)

# [1] - file_path, [2] - heatmap title
arg <- commandArgs(trailingOnly = TRUE)
all_names <- readRDS("all_narna.Rds")
text <- readLines(arg[1])
text <- text[-c(grep("^>", text) %>% max(), length(text))]

names <- sub("^>", "", grep("^>", text, value = TRUE)) # seq names
names <- data.frame("seq_name" = names) %>%
  left_join(all_names, by = "seq_name") %>%
  dplyr::pull(new_name)

data_lines <- text[!grepl("^>", text) & text != ""] # dist values
split_lines <- lapply(data_lines, function(line) as.numeric(strsplit(trimws(line), "\\s+")[[1]]))
dist_values <- unlist(split_lines)
n <- length(split_lines) + 1

dist_matrix <- matrix(0, nrow = n, ncol = n, byrow = TRUE)
idx <- which(lower.tri(dist_matrix), arr.ind = TRUE)
idx_ordered <- idx[order(idx[, 1], idx[, 2]), ]

for (i in seq_len(nrow(idx_ordered))) {
  dist_matrix[idx_ordered[i, "row"], idx_ordered[i, "col"]] <- dist_values[i]
}

dist_matrix <- dist_matrix + t(dist_matrix)
rownames(dist_matrix) <- colnames(dist_matrix) <- names

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

pdf(paste0("hc_and_heatmap/", arg[2], "_heatmap_and_hc.pdf"), width = 12, height = 8)
heatmap(dist_matrix, symm = TRUE, main = paste0("Distance Matrix Heatmap for ", arg[2]), cexRow = 0.5, cexCol = 0.5)
plot(dend, main = paste0("Hierarchical Clustering for ", arg[2]))
dev.off()