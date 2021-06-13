#!/bin/bash
#SBATCH -t 3-00:00:00
#SBATCH -c 1
#SBATCH --mem=80GB

#Now that we have the unsolapped intervals we want the per sample depth corresponding to that intervals.
#With this script we want to run a per species intersect with the samples_bed and the partitioned_bed,
#conserving the intervals of the partitioned_bed and the depth data of each sample.

sp=($(echo $1))

#Copy the intervals file with the name we want to our per species sample coverage bed file.
cp $LUSTRE/${sp}_partitioned.bed $LUSTRE/${sp}_partitioned_samples_coverage.bed

#Create a variable with the sample names
samples=($(ls $STORE2/lynx_genome/lynx_data/CatRef_bams/*cat_ref_sorted_rg_rmdup_sorted_indelrealigner.bam \
| rev | cut -d '/' -f 1 | rev | cut -d '_' -f 1-4))

#With bedtools intersect add the coverage info to the partitioned_bed
for j in ${samples[@]}
  do
  echo "Bedtools intersect coverage for ${j}"
  bedtools intersect -wa -wb \
        -a $LUSTRE/${sp}_partitioned.bed \
        -b $LUSTRE/samples_bed/${j}.coverage.bed | cut -f7 \
        > tmp
  paste $LUSTRE/${sp}_partitioned_samples_coverage.bed tmp > tmp2
  mv tmp2 $LUSTRE/${sp}_partitioned_samples_coverage.bed
  done
  echo "Done bed coverage for ${sp}"
