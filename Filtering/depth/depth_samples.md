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

Now that I have access to CESGA, I will repeat the previous step for the rest of lynx species. The first thing I have to do is copy Felis_catus.100x100kbp.masked.genome.bed to CESGA (from my laptop)

```
scp Felis_catus.100x100kbp.masked.genome.bed csbiellf@ft2.cesga.es:/home/csic/bie/llf
```
Then I will have to run the previous code in an .sh format (See samtools_depth.sh)
```
scp /Users/lorenalorenzo/github/selection_scan_lynx/Filtering/depth/executables/samtools_depth.sh csbiellf@ft2.cesga.es:/home/csic/bie/llf
```
