library(SummarizedExperiment)
library(tidyverse)
library(ConsensusClusterPlus)
library(DESeq2)

salmon_output <- readRDS('all_samples_gene.SummarizedExperiment.rds')
tpms <- assay(salmon_output)
filtered_tpms <- tpms %>% filter((rowSums(tpms==0)/ncol(tpms))<0.3)
scaled_tpms <- log2(filtered_tpms+1) %>% as.matrix()#%>% scale()

top_sd_genes <- rowSds(scaled_tpms) %>% sort(decreasing = T) %>% names()
cutoff <- 2000
top_scaled_tpms <- scaled_tpms[top_sd_genes[1:cutoff],]
top_tpms <- filtered_tpms[top_sd_genes[1:cutoff],]

res1 <- ConsensusClusterPlus(d=as.matrix(top_tpms), maxK = 10, reps=100, pItem=0.8, clusterAlg="hc", title="top_clusters", plot='png')
res2 <- ConsensusClusterPlus(d=as.matrix(top_scaled_tpms), maxK = 10, reps=100, pItem=0.8, clusterAlg="hc", title="sacled_clusters", plot='png')
res3 <- ConsensusClusterPlus(d=as.matrix(filtered_tpms), maxK = 10, reps=100, pItem=0.8, clusterAlg="hc", title="filtered_clusters", plot='png')


km <- kmeans(t(as.matrix(top_tpms)), centers = 4, iter.max = 1000, nstart = 25)
km_scale <- kmeans(t(scaled_tpms[top_sd_genes[1:cutoff],]), centers = 4, iter.max = 1000, nstart = 50)

data.frame(km$cluster, km_scale$cluster) %>% table()
data.frame(res3[[4]]$consensusClass, km_scale$cluster) %>% table()
data.frame(res3[[4]]$consensusClass, km$cluster) %>% table()
klasses <- data.frame(km$cluster, km_scale$cluster, res3[[4]]$consensusClass) %>% set_names(c('top', 'scaled_top', 'cons'))

design <- data.frame(km_scale$cluster) %>% rownames_to_column("Sample") %>% rename("Cluster"='km_scale.cluster') %>% mutate(Cluster=as.factor(Cluster))

#PCA
pca <- prcomp(t(scaled_tpms), scale=T)
pca_data <- as.data.frame(pca$x[, 1:3]) %>% rownames_to_column("Sample") %>%
  left_join(design, by='Sample')
barplot_data <- data.frame(Component = paste0('PC', 1:length(pca$sdev)),
                           Standart_Deviation = pca$sdev)
ggplot(barplot_data[1:9,], aes(x = Component, y = Standart_Deviation)) +
  geom_bar(stat = 'identity', fill = 'lightblue') +
  labs(x = "Principal Component", y = "Standart Deviation") + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
filtered_var <- paste0("PC", 1:length(pca$sdev), ' ', round(pca$sdev / sum(pca$sdev)))
ggplot(pca_data, aes(x = PC1, y = PC2, color = as.factor(Condition))) +
  geom_point(pch = 19, size = 3) +
  labs(x = filtered_var[1], y = filtered_var[2]) +
  theme_minimal() +
  ggtitle("PCA plot") +
  geom_text(aes(label = Sample), vjust = -0.5, hjust = 0.5) + 
  labs(colour = "Cluster")

#DESeq
dds <- filtered_tpms %>% mutate(across(where(is.numeric), as.integer)) %>% as.matrix() %>% DESeqDataSetFromMatrix(., design, design= ~Cluster)
dds <- DESeq(dds)
res <- results(dds) %>% as.data.frame()
top_rank <- res %>% mutate(rank=sign(log2FoldChange)*log10(pvalue)) %>% arrange(desc(rank)) %>% filter(abs(rank)>10)
