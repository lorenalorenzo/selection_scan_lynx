# Purpose: This script is used to split the vcf files by chromosome and species for lassi input
#name variables
species=(lc ll lp lr)
INPUT=$STORE/saltilassi/data
OUTPUT=$STORE/saltilassi/input

for sp in ${species[@]}
do
  echo "$sp"
  #make a sp_ind.txt (pop file) needed for lassip
    sed "s/$/\t${sp}/" <(grep "#CHR" $INPUT/${sp}_goodsamples_filtered_cat_ref.vcf| tr "\t" "\n" | grep ${sp})  > $OUTPUT/${sp}_ind.txt
    sbatch lassi.sh $sp
  CHR=($(cat $STORE2/reference_genomes/Felis_catus_Ref/Felis_catus.Felis_catus_9.0.dna.toplevel.fa.fai | cut -f 1 | uniq))
  for chr in ${CHR[@]:0:18}
    do
      echo "$chr"
      grep -E "^(#|${chr})" \
      $INPUT/${sp}_goodsamples_cat_ref.filter8.vcf \
      > $OUTPUT/${chr}_${sp}_goodsamples_filtered_cat_ref.vcf
      echo "Done $chr for $sp vcf"
    done
done