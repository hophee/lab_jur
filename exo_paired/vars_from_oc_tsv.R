library(tidyverse)
library(janitor)

get_vars_txt <- function(file_path, oc=F, mis=F, file_name){
  data <- read_tsv(file_path, skip=6) %>% janitor::clean_names() %>% select(where(~ !all(is.na(.))))
  if(oc==F){
    if(mis==F){
      paste0(data$gene, ':', data$c_dna_change, "    #", data$vcf_ref_allele, '>',data$vcf_alt_allele) %>% write(file=file_name)
    } else{
      paste0(data$chrom_1, ":", data$position, ' ', data$ref_base, '-', data$alt_base) %>% write(file=file_name)
    }
  } else{
      paste(data$chrom_1, data$position, data$ref_base, data$alt_base) %>%  write(file=file_name)
  }
}
get_vars_txt(file_path="../missense_germ_hard_path.tsv", file_name = 'missens_germ_vars.txt')
get_vars_txt(file_path="../tumor_mis.tsv", file_name = 'missens_tumor_vars.txt')
get_vars_txt(file_path="../lof_germline.tsv", file_name = 'lof_germ_vars.txt')
get_vars_txt(file_path="../hard_tumor_lof.tsv", file_name = 'lof_tumor_vars.txt')