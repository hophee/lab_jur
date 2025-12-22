library(Biostrings)
library(dplyr)
library(readr)

genome_table <- read_tsv('test_data/ECZV.tsv', skip=2)
#genome_table %>% mutate(len=Stop-Start+1) %>% filter(len<=2000 & len >= 1000) %>% arrange(desc(len))
genome <- readDNAStringSet("test_data/ECZV.fna")[[1]]

g_idx <- which(genome_table$Gene == 'tktA')[1]
start <- genome_table$Start[g_idx]
end <- genome_table$Stop[g_idx]

if((end-start+1)<500){
  target <- genome[start:end]
}else{
  mid <- floor((end-start)/2)+start
  side <- (((end-start)/10 + 500)/2) %>% floor()
  target <- genome[(mid-side+1):(mid+side)]
}
