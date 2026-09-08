#!/usr/bin/env bash
# WES pipeline
# Tools: FastQC, HISAT2, samtools, GATK 4.1.9.0, ANNOVAR
set -euo pipefail

# EDIT THESE PATHS / NAMES
SAMPLE="X"
FAMILY="XFamily"
THREADS=24

GATK_JAR=".../gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar"
HISAT2_INDEX=".../hg19/genome"
REF=".../hg19/hg19.fa"

KNOWN_DBSNP=".../Variant_Known_Database_dbsnp138_hg19/dbsnp_138.hg19.vcf"
HAPMAP=".../hapmap_3.3.hg19.vcf"
OMNI=".../1000G_omni2.5.hg19.vcf"
KG_SNP=".../1000G_phase1.snps.high_confidence.hg19.vcf"
CEU=".../CEUTrio.HiSeq.WGS.b37.bestPractices.hg19.vcf"
MILLS=".../Mills_and_1000G_gold_standard.indels.hg19.vcf"
KG_INDEL=".../1000G_phase1.indels.hg19.vcf"

ANNOVAR_PL=".../annovar/table_annovar.pl"
ANNOVAR_XREF=".../annovar/example/gene_fullxref.txt"
ANNOVAR_DB=".../AnnovarDBs"

# FASTQ inputs
FWD_PAIRED="ForwardPairs.fastq"
REV_PAIRED="ReversePairs.fastq"
FWD_UNPAIRED="ForwardUnpairs.fastq"
REV_UNPAIRED="ReverseUnpairs.fastq"

# Extra GVCFs for family joint calling (edit list)
FAMILY_GVCFS=(
  "${SAMPLE}.gvcf"
  "${SAMPLE}.gvcf"
  "${SAMPLE}.gvcf"
  "${SAMPLE}.gvcf"
  "${SAMPLE}.gvcf"
)

gatk() {
  java -jar "${GATK_JAR}" "$@"
}

# FASTQC
fastqc -f fastq *.fastq

# ALIGNING (hisat2)
hisat2 -q -p "${THREADS}" --add-chrname \
  -x "${HISAT2_INDEX}" \
  -1 "${FWD_PAIRED}" \
  -2 "${REV_PAIRED}" \
  -U "${FWD_UNPAIRED}" \
  -U "${REV_UNPAIRED}" \
  -S "${SAMPLE}.sam"

# SAM to BAM & filtration
samtools view -b -F 4 --threads "${THREADS}" "${SAMPLE}.sam" > "${SAMPLE}_unsorted.bam"

# Sorting BAM
gatk SortSam \
  -I "${SAMPLE}_unsorted.bam" \
  -O "${SAMPLE}_Sorted.bam" \
  --SORT_ORDER coordinate

# Read grouping (if forgotten)
gatk AddOrReplaceReadGroups \
  -I "${SAMPLE}_Sorted.bam" \
  -O "${SAMPLE}_Sorted_RG.bam" \
  --RGLB Twist-WES \
  --RGPL ILLUMINA \
  --RGPU MD \
  --RGSM "${SAMPLE}"

# MarkDuplicates
gatk MarkDuplicates \
  -I "${SAMPLE}_Sorted_RG.bam" \
  -M "${SAMPLE}.txt" \
  -O "${SAMPLE}_Sorted_RG_MD.bam"

# BaseRecalibrator
# step1
gatk BaseRecalibrator \
  -I "${SAMPLE}_Sorted_RG_MD.bam" \
  -O "${SAMPLE}BQSR.txt" \
  -R "${REF}" \
  --known-sites "${KNOWN_DBSNP}"

# step2
gatk ApplyBQSR \
  -bqsr "${SAMPLE}BQSR.txt" \
  -I "${SAMPLE}_Sorted_RG_MD.bam" \
  -O "${SAMPLE}_Sorted_RG_MD_bqsr.bam"

gatk BaseRecalibrator \
  -I "${SAMPLE}_Sorted_RG_MD_bqsr.bam" \
  -O "${SAMPLE}_table.txt" \
  -R "${REF}" \
  --known-sites "${KNOWN_DBSNP}"

# AnalyzeCovariates
gatk AnalyzeCovariates \
  -before "${SAMPLE}BQSR.txt" \
  -after "${SAMPLE}_table.txt" \
  -plots "${SAMPLE}_plots.pdf" \
  -csv "${SAMPLE}_CSV.csv"

# Joint-Calling methods of multiple DNA samples via GVCF workflow
# HaplotypeCaller
gatk HaplotypeCaller \
  -I "${SAMPLE}_Sorted_RG_MD_bqsr.bam" \
  -O "${SAMPLE}.gvcf" \
  -ERC GVCF \
  -R "${REF}" \
  --sample-name "${SAMPLE}"

# combining
COMBINE_ARGS=()
for v in "${FAMILY_GVCFS[@]}"; do
  COMBINE_ARGS+=(-V "${v}")
done

gatk CombineGVCFs \
  -R "${REF}" \
  "${COMBINE_ARGS[@]}" \
  -O "${FAMILY}_combineHC.vcf"

# joining
gatk GenotypeGVCFs \
  -R "${REF}" \
  -V "${FAMILY}_combineHC.vcf" \
  -O "${FAMILY}_joinHC.vcf"

# Counting Variants
gatk CountVariants \
  -V "${FAMILY}_joinHC.vcf"

# VQSR for SNPs
# step1
gatk VariantRecalibrator \
  -O "${FAMILY}_SNP_Recal.vcf" \
  --resource:hapmap,known=false,training=true,truth=true,prior=15.0 "${HAPMAP}" \
  --resource:omni,known=false,training=true,truth=false,prior=12.0 "${OMNI}" \
  --resource:1000G,known=false,training=true,truth=false,prior=10.0 "${KG_SNP}" \
  --resource:dbSNP_138,known=true,training=false,truth=false,prior=2.0 "${KNOWN_DBSNP}" \
  --resource:CEU,known=true,training=false,truth=false,prior=2.0 "${CEU}" \
  --tranches-file "TranchesSNP.txt" \
  -an QD -an FS -an SOR -an ReadPosRankSum \
  -V "${FAMILY}_joinHC.vcf" \
  --max-gaussians 8 \
  -R "${REF}" \
  --rscript-file "RscriptSNP.R" \
  --mode SNP \
  -tranche 90.0 -tranche 91.0 -tranche 92.0 -tranche 93.0 -tranche 94.0 \
  -tranche 95.0 -tranche 96.0 -tranche 97.0 -tranche 98.0 -tranche 99.0 \
  -tranche 99.9 -tranche 100.0

# step2
gatk ApplyVQSR \
  --recal-file "${FAMILY}_SNP_Recal.vcf" \
  -V "${FAMILY}_joinHC.vcf" \
  -O "${FAMILY}_HC_SNP_VQSR_done.vcf" \
  -mode SNP \
  -R "${REF}" \
  --tranches-file "TranchesSNP.txt" \
  -ts-filter-level 99.9

# VQSR for Indels
# step1
gatk VariantRecalibrator \
  -O "${FAMILY}_HC_Indels_Recal.vcf" \
  --resource:Mills,known=false,training=true,truth=true,prior=12.0 "${MILLS}" \
  --resource:dbSNP_138,known=true,training=false,truth=false,prior=2.0 "${KNOWN_DBSNP}" \
  --resource:CEU,known=true,training=false,truth=false,prior=2.0 "${CEU}" \
  --resource:1000G,known=false,training=true,truth=false,prior=10.0 "${KG_INDEL}" \
  --tranches-file "${FAMILY}.Indel.Tranche" \
  -an QD -an FS -an SOR -an ReadPosRankSum \
  -V "${FAMILY}_HC_SNP_VQSR_done.vcf" \
  --max-gaussians 8 \
  -R "${REF}" \
  --rscript-file "${FAMILY}-rscript.Indels.R" \
  --mode INDEL \
  -tranche 90.0 -tranche 91.0 -tranche 92.0 -tranche 93.0 -tranche 94.0 \
  -tranche 95.0 -tranche 96.0 -tranche 97.0 -tranche 98.0 -tranche 99.0 \
  -tranche 99.9 -tranche 100.0

# step2
gatk ApplyVQSR \
  --recal-file "${FAMILY}_HC_Indels_Recal.vcf" \
  -V "${FAMILY}_HC_SNP_VQSR_done.vcf" \
  -O "${FAMILY}_HC_SNPs_INDELS_VQSR_done.vcf" \
  -mode INDEL \
  -R "${REF}" \
  --tranches-file "${FAMILY}.Indel.Tranche" \
  -ts-filter-level 99.9

# ANNOTATION
perl "${ANNOVAR_PL}" \
  --xreffile "${ANNOVAR_XREF}" \
  --polish \
  --thread 24 \
  --dot2underline \
  --nastring . \
  --vcfinput \
  --otherinfo \
  --remove \
  --buildver hg19 \
  --outfile "${FAMILY}_annotated" \
  --protocol refGeneWithVer,ensGene,cytoBand,genomicSuperDups,clinvar_20210123,1000g2015aug_all,gnomad211_genome,gnomad211_exome,exac03,kaviar_20150923,esp6500siv2_all,gme,hrcr1,iranome,Coding_CosmicV92,Noncoding_CosmicV92,cg46,cg69,icgc28,intervar_20180118,dann,cadd13,CADD16,dbnsfp41a,avsnp144,avsnp150,fathmm \
  --operation g,g,r,r,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f \
  "${FAMILY}_HC_SNPs_INDELS_VQSR_done.vcf" \
  "${ANNOVAR_DB}"
