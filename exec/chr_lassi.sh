#!/bin/bash
#SBATCH --time=3-00:00:00
#SBATCH --mem=10GB

#name variables
sp=($(echo $1))
chr=($(echo $2))
INPUT=$STORE/saltilassi/${chr}_${sp}_goodsamples_filtered_polarized_variants_header_cat_ref.vcf
#MAP=$STORE2/reference_genomes/Felis_catus_Ref/Felis_catus.Felis_catus_9.0.dna.toplevel.fa.fai
#POP=$LUSTRE/selection_scan/${sp}_ind.txt
#OUT=$LUSTRE/selection_scan/saltiLASSI/${chr}_${sp}_salti

#load dependencies
module load lassip

#run lassi
echo "running salti-lassi in $chr"    
lassip \
          --vcf INPUT \
          --unphased \
          --calc-spec \
          --hapstats \
          --salti \
          --map MAP \
          --winsize 101 \
          --winstep 50 \
          --pop POP \
          --out OUT

