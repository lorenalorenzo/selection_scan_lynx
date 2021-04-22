---
Title: Calculating mean depth per sample
Author: Lorena Lorenzo Fernández
Date: 02 April, 2021
---
The whole genome is too big to analyze by region depth, that's why we only need a representative subset of data. For doing this, I am going to subsample BAM files, previously generating a BED file with 100 random positions of 100000 bp length. Thus, depth calculations will be done only considering this random subset.

To generate the random regions file I will use BEDtools random. I'll then remove the low-mappability and repetitive regions from the file using BEDtools subtract

```
# Create a Genome region file for Bedtools:
# A file with the list of chromosomes as col1 and their length as col2, tab separated
# Basically the first two columns of a FAI file:
# Using an fai index file in conjunction with a FASTA/FASTQ file containing reference sequences enables efficient access to arbitrary regions within those reference sequences. The index file typically has the same filename as the corresponding FASTA/FASTQ file, with .fai appended.

cut -f1,2 /GRUPOS/grupolince/reference_genomes/felis_catus_genome/Felis_catus.Felis_catus_9.0.dna.toplevel.fa.fai > /home/llorenzo/depth_filter/Felis_catus_9.0.dna.toplevel.genome

# Bedtools random to generate file of 100 random segments of 100000 bp
# Output a BED file:

bedtools random -l 100000 -n 100 -g /home/llorenzo/depth_filter/Felis_catus_9.0.dna.toplevel.genome | sort > /home/llorenzo/depth_filter/Felis_catus.100x100kbp.genome.bed

# Using bedtools subtract I can remove low-mappability and repetitive regions:

bedtools subtract -a /home/llorenzo/depth_filter/Felis_catus.100x100kbp.genome.bed -b /GRUPOS/grupolince/reference_genomes/felis_catus_genome/Masked_Regions.bed > /home/llorenzo/depth_filter/Felis_catus.100x100kbp.masked.genome.bed
```

Then, I have to calculate region coverage per sample with samtools depth (loop to do this in every sample). Depth at all positions will be calculated (-a) within the regions randomly selected before (-b). What I am going to use is bams of lp until I get access to CESGA (in process) and to the lynx table (Dropbox)

This will be run in a loop for each lp bam sample.

```
# Loop of Samtools depth calculations for each sample_bam
for i in *lp*.bam
  do
  echo "Calculating depth for $i"
  samtools depth -a -b /home/llorenzo/depth_filter/Felis_catus.100x100kbp.masked.genome.bed \
  "$i" \
  > /home/llorenzo/depth_filter/"$i".100x100kbp.masked.depth
done
```

With this, I will analyze coverage distribution (average, stdv, maxDepth, minDepth) and summarize information in a table. See depth_script.R

```
scp llorenzo@genomics-a.ebd.csic.es:/home/llorenzo/depth_filter/*masked.depth .
```


## Once I had access to CESGA, REST OF THE SPECIES

Now that I have access to CESGA, I will repeat the previous step for the rest of lynx species. The first thing I have to do is copy Felis_catus.100x100kbp.masked.genome.bed to CESGA (from my laptop)

```
scp Felis_catus.100x100kbp.masked.genome.bed csbiellf@ft2.cesga.es:/home/csic/bie/llf
```
Then I will have to run the previous code in an .sh format (See samtools_depth.sh). Upload the sh in the CESGA server:
```
scp /Users/lorenalorenzo/github/selection_scan_lynx/Filtering/depth/executables/samtools_depth.sh csbiellf@ft2.cesga.es:/home/csic/bie/llf
```
and run it:
```
sbatch samtool_depth.sh
```
In this way, I estimated 16h to do this work (obviously not enough), so only lc samples are calculated (except one: c_lc_zz_0003)
As doing the loop in every sample will be high time consuming (more than 3 days), I will loop the sbatch for every sample, in other words, send 81 individual works (as we have 81 samples) to the server.
Copy the new sh from the laptop to the server:
```
scp /Users/lorenalorenzo/github/selection_scan_lynx/Filtering/depth/executables/samtools_depth_per_sample.sh csbiellf@ft2.cesga.es:/home/csic/bie/llf
```
Nevertheless, as lp samples are already calculated and lc as well (except the one mentioned above), I do not calculate on this samples to make the process faster and do not duplicate results. So, from 81 we pass to 81-11lp-18lc=52 samples(lr,ll and one lc). The lc sample will be send alone as it is difficult for me to include in the list of $sample. In adittion, "c_ll_ca_0249" and "c_ll_ca_0253" are bad samples which are not counted in the total samples (so, 81 samples without the bad ones, remaining 52 samples).

So, for the rest:

```
sample=($(ls $STORE2/lynx_genome/lynx_data/CatRef_bams/*cat_ref_sorted_rg_rmdup_sorted_indelrealigner.bam \
| grep -v "lp" | grep -v "lc" | rev | cut -d '/' -f 1 | rev | cut -d '_' -f 1-4))

for i in ${sample[@]}
do
echo "Calculating depth for $i"
sbatch samtools_depth_per_sample.sh $i
done
```
And for c_lc_zz_0003:
```
lcsample=($(ls $STORE2/lynx_genome/lynx_data/CatRef_bams/c_lc_zz_0003*cat_ref_sorted_rg_rmdup_sorted_indelrealigner.bam \
| rev | cut -d '/' -f 1 | rev | cut -d '_' -f 1-4))

for i in ${lcsample}
do
echo "Calculating depth for $i"
sbatch samtools_depth_per_sample.sh $lcsample
done
```
Download results in my laptop for R analyses:
```
scp -r csbiellf@ft2.cesga.es:/mnt/lustre/scratch/home/csic/bie/llf .
```
