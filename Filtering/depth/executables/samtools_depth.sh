module load Samtools

# Loop of Samtools depth calculations for each sample_bam (only for lp)
for i in *cat_ref_sorted_rg_rmdup_sorted_indelrealigner.bam
  do
  echo "Calculating depth for $i"
  samtools depth -a -b /home/csic/bie/llf/Felis_catus.100x100kbp.masked.genome.bed \
  "$i" \
  > /home/csic/bie/llf/"$i".100x100kbp.masked.depth
done

wd of BAM data: /mnt/netapp1/Store_csebdjgl/lynx_genome/lynx_data/CatRef_bams
-sbatch!!
