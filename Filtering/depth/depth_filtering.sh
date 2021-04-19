---
Title: Applying depth filter to VCF
Author: Lorena Lorenzo Fernández
Date: 19 April 2021
---

# Now to assign "not-callable" status to <MinDepth and >MaxDepth regions

 for sample in $(ls *.coverage.bed | cut -d'.' -f1)
  do
   echo ${sample}
   max=($(grep "${sample}" 8samples_depth.csv | cut -d',' -f7))
   awk -v max="${max}" '{FS="\t"; OFS="\t"; if ($4 < 6 || $4 > max) print $1, $2, $3, "not-callable"; else print $1, $2, $3, "callable";}' ${sample}.coverage.bed \
   > ${sample}.coverage.defined.bed
 done

 # I can now remove the not-callable regions and also the repetitive/low-mappability ones,
 # to obtain the callable regions of each sample:

 for sample in $(ls *.coverage.defined.bed | cut -d'.' -f1)
  do
   echo ${sample}
   grep -v "not-callable" ${sample}.coverage.defined.bed |
   bedtools subtract -a - -b /GRUPOS/grupolince/reference_genomes/felis_catus_genome/Masked_Regions.bed |
   bedtools merge -i - \
   > ${sample}.callable.bed
 done

 for sample in $(ls *.coverage.defined.bed | cut -d'.' -f1)
  do
   bedtools merge -i ${sample}.callable.bed > ${sample}.callable.collapsed.bed
 done

 # To check length of sequence in BED
 # awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}' c_lc_zz_0001.callable.collapsed.bed

 # Intersect all of the BED files
 multiIntersectBed -i c_lc_zz_0001.callable.collapsed.bed \
 c_lc_zz_0003.callable.collapsed.bed \
 c_ll_vl_0112.callable.collapsed.bed \
 c_ll_ya_0146.callable.collapsed.bed \
 c_lp_sm_0138.callable.collapsed.bed \
 c_lp_sm_0140.callable.collapsed.bed \
 c_lr_nm_0006.callable.collapsed.bed \
 c_lr_zz_0001.callable.collapsed.bed \
 > 8samples.callable.collapsed.intersect.bed

 # Only regions overlapping and callable in all samples
 grep "1,2,3,4,5,6,7,8" 8samples.callable.collapsed.intersect.bed > 8samples.callable.collapsed.intersect.all.bed

 # Get VCF of 8samples from allsamples VCF
 samplesARRAY=($(ls *.coverage.defined.bed | cut -d'.' -f1))

 /opt/gatk-4.1.0.0/gatk SelectVariants \
   -R /GRUPOS/grupolince/reference_genomes/felis_catus_genome/Felis_catus.Felis_catus_9.0.dna.toplevel.fa \
   -V ~/CatRef_vcfs/allsamples_cat_ref.filter5.vcf \
   $(for j in ${samplesARRAY[@]}; do echo "-sn ${j}";done) \
   -O 8samples_cat_ref.filter5.vcf


 # Now we can filter the VCF
 bedtools subtract -a 8samples_cat_ref.filter5.vcf \
   -b 8samples.callable.collapsed.intersect.all.bed \
   -header > 8samples_cat_ref.filter-callable.vcf
