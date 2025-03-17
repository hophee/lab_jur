#!/usr/bin/bash

eval "$(conda shell.bash hook)"

conda activate bwabwa
bwa mem -t 16 Homo_sapiens_assembly38.fasta exome_normal/[sample_name].fastq.gz exome_normal/[sample_name].fastq.gz | samtools sort --threads 3 -m 1G -o normal.bam
samtools index normal.bam
conda deactivate

~/gatk4/gatk MarkDuplicates \
-I normal.bam \
-O normal_dupmark.bam \
-M metrics.txt

~/gatk4/gatk AddOrReplaceReadGroups \
-I normal_dupmark.bam \
-O normal_with_rg.bam \
-LB lib1 \
-PL ILLUMINA \
-PU HXXXLXX \
-SM sample1 \
-CREATE_INDEX true

~/gatk4/gatk BaseRecalibrator \
-R Homo_sapiens_assembly38.fasta \
-I normal_with_rg.bam \
--known-sites vars/Homo_sapiens_assembly38.dbsnp138.vcf \
--known-sites vars/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
-O recalibration.table \
-L bedi/NEXome_Plus_Panel_v1.probe.hg38.bed

~/gatk4/gatk ApplyBQSR \
-R Homo_sapiens_assembly38.fasta \
-I normal_with_rg.bam \
--bqsr-recal-file recalibration.table \
-O recalibrated.bam

~/gatk4/gatk HaplotypeCaller \
-R Homo_sapiens_assembly38.fasta \
-I recalibrated.bam \
-O raw_variants.vcf \
-L bedi/NEXome_Plus_Panel_v1.probe.hg38.bed \
-D vars/Homo_sapiens_assembly38.dbsnp138.vcf \
-G StandardAnnotation

conda activate freebay
bcftools filter -i 'QUAL >= 50 & INFO/DP >= 35 & QD > 2 & MQ > 40.0' raw_variants.vcf -o filtered.vcf
conda deactivate