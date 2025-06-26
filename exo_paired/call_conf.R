library(vcfR)
library(tidyverse)

strelka <- read.vcfR("strelka_merged.vcf.gz")
mutect <- read.vcfR('tumor_sample_vs_normal_sample.mutect2.filtered_snpEff_VEP.ann.vcf.gz')
all <- read.vcfR('all_merge_filtered.vcf.gz')

in_mut <- all@fix[,2] %in% mutect@fix[,2]
in_strelka <- all@fix[,2] %in% strelka@fix[,2]
in_two <- as.logical(in_mut*in_strelka)
all@fix %>% as.data.frame() %>% filter(in_two) %>% pull(POS)
sum(in_mut, in_strelka)

all@fix %>% as.data.frame() %>% pull(POS) %>%  table() %>% sort(decreasing = T)
