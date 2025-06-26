library(vcfR)
library(tidyverse)


hap <- read.vcfR('<filtered_vars_from_nf/rnavar>', verbose = F) #rna_seq_filtered
gt_hap <- extract_gt_tidy(hap)

TDP_hap <- str_split_fixed(gt_hap$gt_AD, ',', 3)[,2:3] %>% as.data.frame() %>% mutate(across(where(is.character), ~ as.numeric(.x) %>% replace_na(0))) %>% rowSums()
CDP_hap <- str_split_fixed(gt_hap$gt_AD, ',', 2)[,1] %>%  as.numeric()
dp <- gt_hap %>% pull(gt_DP)
TVAF_hap <- TDP_hap / dp
CVAF_hap <- CDP_hap / dp

hap@meta <- c(hap@meta , '##INFO=<ID=TVAF,Number=1,Type=Float,Description="Allelic fraction of alternative allele in tumor">', 
                 '##INFO=<ID=CVAF,Number=1,Type=Float,Description="Allelic fraction of alternative allele in control">',
                 '##INFO=<ID=TDP,Number=1,Type=Integer,Description="Tumor sample depth">',
                 '##INFO=<ID=CDP,Number=1,Type=Integer,Description="Control sample depth">')
hap@fix[,8] <- paste0('TVAF=', TVAF_hap, ';', 'CVAF=', CVAF_hap, ';',
                         'TDP=',TDP_hap, ';', 'CDP=', CDP_hap ,';',
                         hap@fix[,8])

write.vcf(hap, file = 'rnaseq.vcf.gz')