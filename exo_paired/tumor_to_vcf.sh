#!/usr/bin/bash

eval "$(conda shell.bash hook)"

conda activate bwabwa
bwa mem -t 16 Homo_sapiens_assembly38.fasta exome_tumor/[tumor_sample_name].fastq.gz exome_tumor/[tumor_sample_name].fastq.gz | samtools sort --threads 4 -m 1G -o tumor.bam
samtools index tumor.bam
conda deactivate

~/gatk4/gatk --java-options "-Xmx10g" MarkDuplicates \
-I tumor.bam \
-O tumor_dupmark.bam \
-M metrics_tumor.txt \
-TMP_DIR ~/tempest

~/gatk4/gatk AddOrReplaceReadGroups \
-I tumor_dupmark.bam \
-O tumor_with_rg.bam \
-LB lib1 \
-PL ILLUMINA \
-PU HXXXLXX \
-SM sample2 \
-CREATE_INDEX true

~/gatk4/gatk --java-options "-Xmx16G" Mutect2 \
-R  Homo_sapiens_assembly38.fasta \
-I tumor_with_rg.bam \
-I normal_with_rg.bam \
-normal sample1 \
--germline-resource vars/af-only-gnomad.hg38.vcf.gz \
--panel-of-normals vars/1000g_pon.hg38.vcf.gz \
-O somatic_variants.vcf.gz \
-L bedi/NEXome_Plus_Panel_v1.probe.hg38.bed \
--native-pair-hmm-threads 16

bcftools filter -i 'TLOD >= 10 & NLOD < 10 & INFO/DP > 200 & POPAF > 2 & PON!=1' somatic_variants.vcf.gz -o filtered_tumor.vcf.gz
