#!/bin/bash
#SBATCH -t 6:00:00
#SBATCH -c 1
#SBATCH --mail-user=lorena.lorenzo.fdez@gmail.com

module load bedtools

sample=($(echo $1))

bedtools genomecov \
-ibam $STORE2/lynx_genome/lynx_data/CatRef_bams/${sample}_cat_ref_sorted_rg_rmdup_sorted_indelrealigner.bam -bga \
> $LUSTRE/"${sample}".coverage.bed
