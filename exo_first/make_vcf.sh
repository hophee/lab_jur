#!/usr/bin/bash

eval "$(conda shell.bash hook)"

conda activate bwabwa
bwa mem -t 16 Homo_sapiens_assembly38.fasta [exome_path1] [exome_path2] | samtools sort --threads 3 -m 1G -o ZD210122.bam
samtools index ZD210122.bam
conda deactivate

conda activate gataka

~/gatk4/gatk MarkDuplicates \
I=ZD210122.bam \
O=ZD210122_dupmark.bam \
M=metrics.txt

~/gatk4/gatk AddOrReplaceReadGroups \
-I ZD210122_dupmark.bam \
-O ZD210122_with_rg.bam \
-LB lib1 \
-PL ILLUMINA \
-PU HXXXLXX \
-SM sample1 \
-CREATE_INDEX true

~/gatk4/gatk BaseRecalibrator \
-R Homo_sapiens_assembly38.fasta \
-I ZD210122_with_rg.bam \
--known-sites vars/Homo_sapiens_assembly38.dbsnp138.vcf \
--known-sites vars/Mills_and_1000G_gold_standard.indels.hg38.vcf \
-O recalibration.table \
-L bedi/S07604514_Covered.bed

~/gatk4/gatk ApplyBQSR \
-R Homo_sapiens_assembly38.fasta \
-I ZD210122_with_rg.bam \
--bqsr-recal-file recalibration.table \
-O recalibrated.bam


~/gatk4/gatk HaplotypeCaller \
-R Homo_sapiens_assembly38.fasta \
-I recalibrated.bam \
-O raw_variants.vcf \
-L bedi/S07604514_Covered.bed \
-D vars/Homo_sapiens_assembly38.dbsnp138.vcf \
-G StandardAnnotation

bcftools filter -i 'QUAL >= 50 & INFO/DP >= 35 & QD > 2 & MQ > 40.0' raw_variants.vcf -o filtered.vcf
#полученные на этом этапе vcf загружаются в wannovar и vep
#обрабатываются дальеш скриптом filter.R

#дальше идёт неиспользуемая часть решения
#стоит удалить эту часть перед запуском
~/gatk4/gatk VariantRecalibrator \
-R Homo_sapiens_assembly38.fasta \
-V raw_variants.vcf \
--resource:hapmap,known=false,training=true,truth=true,prior=15.0 known_sites/hapmap_3.3.hg38.vcf.gz \
--resource:omni,known=false,training=true,truth=false,prior=12.0 known_sites/1000G_omni2.5.hg38.vcf.gz \
--resource:1000G,known=false,training=true,truth=false,prior=10.0 known_sites/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
--resource:dbsnp,known=true,training=false,truth=false,prior=2.0 vars/Homo_sapiens_assembly38.dbsnp138.vcf \
-an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
--mode SNP \
--output variants_SNP.recal \
--tranches-file variants_SNP.tranches

~/gatk4/gatk VariantRecalibrator \
-R Homo_sapiens_assembly38.fasta \
-V raw_variants.vcf \
--resource:mills,known=true,training=true,truth=true,prior=12.0 vars/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
--resource:dbsnp,known=true,training=false,truth=false,prior=2.0 vars/Homo_sapiens_assembly38.dbsnp138.vcf \
-an QD -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
--mode INDEL \
--output variants_INDEL.recal \
--tranches-file variants_INDEL.tranches

# 4. Применение VQSR для обоих типов вариантов
~/gatk4/gatk ApplyVQSR \
-R Homo_sapiens_assembly38.fasta \
-V raw_variants.vcf \
--recal-file variants_SNP.recal \
--tranches-file variants_SNP.tranches \
--truth-sensitivity-filter-level 99.0 \
--mode SNP \
-O recalibrated_SNP.vcf

~/gatk4/gatk ApplyVQSR \
-R Homo_sapiens_assembly38.fasta \
-V recalibrated_SNP.vcf \
--recal-file variants_INDEL.recal \
--tranches-file variants_INDEL.tranches \
--truth-sensitivity-filter-level 99.0 \
--mode INDEL \
-O recalibrated_variants.vcf  # Финальный файл
conda deactivate

conda activate snef
# 5. Фильтрация PASS-вариантов
SnpSift filter -f recalibrated_variants.vcf "FILTER = 'PASS'" > recalibrated_pass.vcf


# Аннотация
snpEff -Xmx8G -v GRCh38.99 -stats snpEff_stats.html recalibrated_pass.vcf > annotated_variants.vcf 2> snpEff.log

# Добавление ClinVar
SnpSift annotate -info CLNSIG,CLNDN,CLNREVSTAT vars/clinvar.vcf.gz annotated_variants.vcf > annotated_with_clinvar.vcf

# Фильтрация по воздействию
SnpSift filter \
"(ANN[*].IMPACT == 'HIGH') | (ANN[*].IMPACT == 'MODERATE')" \
annotated_with_clinvar.vcf \
> VARIANTS_VQSRed_pass_clinvar_impact.vcf

# Поиск патогенных вариантов
grep -v "##" VARIANTS_VQSRed_pass_clinvar_impact.vcf | grep "pathogenic" > pathogenic_vars.txt

conda deactivate
