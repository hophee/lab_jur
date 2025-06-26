library(vcfR)
library(tidyverse)

snps <- read.vcfR('tumor_sample_vs_normal_sample.strelka.somatic_snvs_snpEff_VEP.ann.vcf.gz', verbose = F)
#vcfR2tidy() for work with df
gt_snps <- extract_gt_tidy(snps, verbose = F)
#REF allele
vcf_ref_U <- paste0(snps@fix[,4], 'U') %>% tolower()
#ALT allele
vcf_alt_U <- paste0(snps@fix[,5], 'U') %>% tolower()

#freq of tier1 alleles
au <- str_split_fixed(gt_snps$gt_AU[gt_snps$Indiv=='TUMOR'], n=2, ',')[,1] %>% as.numeric()
tu <- str_split_fixed(gt_snps$gt_TU[gt_snps$Indiv=='TUMOR'], n=2, ',')[,1] %>% as.numeric()
gu <- str_split_fixed(gt_snps$gt_GU[gt_snps$Indiv=='TUMOR'], n=2, ',')[,1] %>% as.numeric()
cu <- str_split_fixed(gt_snps$gt_CU[gt_snps$Indiv=='TUMOR'], n=2, ',')[,1] %>% as.numeric()
af_tier1_snps <- cbind(au, tu, gu, cu)
rm(au, tu, gu, cu)

#calculate AF. WARNING: required a lot of RAM. Very inefficient code.
gc(verbose=F)
ref_af_count <- diag(af_tier1_snps[1:dim(af_tier1_snps)[1], vcf_ref_U])
gc(verbose=F)
alt_af_count <- diag(af_tier1_snps[1:dim(af_tier1_snps)[1], vcf_alt_U])
gc(verbose=F)

AF_snps <- alt_af_count/(alt_af_count+ref_af_count)


#indels  AF count
indels <- read.vcfR('tumor_sample_vs_normal_sample.strelka.somatic_indels_snpEff_VEP.ann.vcf.gz', verbose = F)
gt_indels  <- extract_gt_tidy(indels, verbose = F)
tar <- str_split_fixed(gt_indels$gt_TAR[gt_indels$Indiv=='TUMOR'], n=2, ',')[,1] %>% as.numeric()
tir <- str_split_fixed(gt_indels$gt_TIR[gt_indels$Indiv=='TUMOR'], n=2, ',')[,1] %>% as.numeric()
AF_indels <- tir/(tar+tir)
rm(tar,tir)

#get control tags
au <- str_split_fixed(gt_snps$gt_AU[gt_snps$Indiv=='NORMAL'], n=2, ',')[,1] %>% as.numeric()
tu <- str_split_fixed(gt_snps$gt_TU[gt_snps$Indiv=='NORMAL'], n=2, ',')[,1] %>% as.numeric()
gu <- str_split_fixed(gt_snps$gt_GU[gt_snps$Indiv=='NORMAL'], n=2, ',')[,1] %>% as.numeric()
cu <- str_split_fixed(gt_snps$gt_CU[gt_snps$Indiv=='NORMAL'], n=2, ',')[,1] %>% as.numeric()
af_tier1_snps_control <- cbind(au, tu, gu, cu)
rm(au, tu, gu, cu)
gc()
control_ref_af_count <- diag(af_tier1_snps_control[1:dim(af_tier1_snps_control)[1], vcf_ref_U])
gc()
control_alt_af_count <- diag(af_tier1_snps_control[1:dim(af_tier1_snps_control)[1], vcf_alt_U])
gc()

control_AF_snps <- control_alt_af_count/(control_alt_af_count+control_ref_af_count)

control_tar <- str_split_fixed(gt_indels$gt_TAR[gt_indels$Indiv=='NORMAL'], n=2, ',')[,1] %>% as.numeric()
control_tir <- str_split_fixed(gt_indels$gt_TIR[gt_indels$Indiv=='NORMAL'], n=2, ',')[,1] %>% as.numeric()
control_AF_indels <- control_tir/(control_tar+control_tir)
rm(control_tir, control_tar)

#get DP
tdp_indels <- gt_indels$gt_DP[gt_indels$Indiv=='TUMOR']
cdp_indels <- gt_indels$gt_DP[gt_indels$Indiv=='NORMAL']
tdp_snps <- gt_snps$gt_DP[gt_snps$Indiv=='TUMOR']
cdp_snps <- gt_snps$gt_DP[gt_snps$Indiv=='NORMAL']


#Write
#snps@meta <- c(snps@meta, '##INFO=<ID=TVAF,Number=1,Type=Float,Description="Allelic fraction of alternative allele in tumor">',
#               '##FORMAT=<ID=AL,Number=1,Type=Integer,Description="Codes of alghoritm  that called somatic call (mutect = 1, strelka = 2)">',
#               '##INFO=<ID=TAL,Number=.,Type=String,Description="Algorithms that called the somatic mutation">')


indels@fix[,8] <- paste0('TVAF=', round(AF_indels,3), ';', 'CVAF=', round(control_AF_indels,3), ';',
                         'TDP=',tdp_indels, ';', 'CDP=', cdp_indels ,';', 
                         indels@fix[,8])
indels@meta <- c(indels@meta, '##INFO=<ID=TVAF,Number=1,Type=Float,Description="Allelic fraction of alternative allele in tumor">', 
                 '##INFO=<ID=CVAF,Number=1,Type=Float,Description="Allelic fraction of alternative allele in control">',
                 '##INFO=<ID=TDP,Number=1,Type=Integer,Description="Tumor sample depth">',
                 '##INFO=<ID=CDP,Number=1,Type=Integer,Description="Control sample depth">')
write.vcf(indels, file='strelka_indels.vcf.gz')

snps@fix[,8] <- paste0('TVAF=', round(AF_snps,3), ';', 'CVAF=', round(control_AF_snps,3), ';',
                       'TDP=',tdp_snps, ';', 'CDP=', cdp_snps ,';',
                       snps@fix[,8])
snps@meta <- c(snps@meta, '##INFO=<ID=TVAF,Number=1,Type=Float,Description="Allelic fraction of alternative allele in tumor">', 
                 '##INFO=<ID=CVAF,Number=1,Type=Float,Description="Allelic fraction of alternative allele in control">',
               '##INFO=<ID=TDP,Number=1,Type=Integer,Description="Tumor sample depth">',
               '##INFO=<ID=CDP,Number=1,Type=Integer,Description="Control sample depth">')
write.vcf(snps, file='strelka_snps.vcf.gz')


#FOR MUTECT GET AF AND DP
mutect <- read.vcfR('tumor_sample_vs_normal_sample.mutect2.filtered_snpEff_VEP.ann.vcf.gz')
gt_mut <- extract_gt_tidy(mutect)
info_mut <- extract_info_tidy(mutect)
TVAF_mut <- gt_mut %>% filter(Indiv=='patient1_tumor_sample') %>% pull(gt_AF) 
CVAF_mut <- gt_mut %>% filter(Indiv=='patient1_normal_sample') %>% pull(gt_AF) 
TDP_mut <- gt_mut %>% filter(Indiv=='patient1_tumor_sample') %>% pull(gt_DP) 
CDP_mut <- gt_mut %>% filter(Indiv=='patient1_normal_sample') %>% pull(gt_DP)
mutect@meta <- c(mutect@meta , '##INFO=<ID=TVAF,Number=1,Type=Float,Description="Allelic fraction of alternative allele in tumor">', 
                 '##INFO=<ID=CVAF,Number=1,Type=Float,Description="Allelic fraction of alternative allele in control">',
                 '##INFO=<ID=TDP,Number=1,Type=Integer,Description="Tumor sample depth">',
                 '##INFO=<ID=CDP,Number=1,Type=Integer,Description="Control sample depth">')
mutect@fix[,8] <- paste0('TVAF=', TVAF_mut, ';', 'CVAF=', CVAF_mut, ';',
                     'TDP=',TDP_mut, ';', 'CDP=', CDP_mut ,';',
                     mutect@fix[,8])
colnames(mutect@gt) <- colnames(snps@gt)
write.vcf(mutect, file = 'mutect.vcf.gz')



