---
Title: Applying depth filter to VCF
Author: Lorena Lorenzo Fernández
Date: 19 April 2021
---
At this point we may have evaluated per sample and species depth distribution (See depth_sample_distribution.R and depth_species_distribution.R). We have several ways to use this info:
    A. Use min. and max. depth generated from distribution to filter each sample BED.
    B. Eliminate samples with low quality of incongruences and then filter min and max depth separately
We are not going to follow the A option because that way we will discard much "good info".

## B. Eliminate samples with low quality of incongruences and then filter min and max depth separately
The first thing we are going to do is discard one sample: c_lc_yu_0007

#### MÁXIMUM depth filter
We decided to analyze depth distribution per species in order to avoid duplications but keep good information that may be isn't present in every sample due to quality issues. So that filter is a trade off.

First, we need to join each sample_depth info per species (remember the sample called c_lr_xx_0011 is a lc) and then analyze distribution in R similarly as done before with each sample.
```
#in CESGA
cd $LUSTRE/samples_subset_depth
cut -f1,2 c_lc_ak_0015.100x100kbp.masked.depth > depth_rows

species=(lc ll lp lr)

for i in ${species[@]}
 do
  echo ${i}
  samples=($(ls $LUSTRE/samples_subset_depth/*masked.depth |rev | cut -d '/' -f 1 | grep "${i}" | rev | cut -d '.' -f 1)))
  for j in ${samples[@]}
   do
    echo ${j}
    cut -f3 ${j}.100x100kbp.masked.depth > ${j}.depthcolumn
  done
  paste depth_rows *${i}*.depthcolumn > ${i}_allsamples.depth
done
```
Now I have to copy the outputs to my laptop in order to analyze in R:
```
scp csbiellf@ft2.cesga.es:/mnt/lustre/scratch/home/csic/bie/llf/samples_subset_depth/\*_allsamples.depth .
```
I had to write a backlash (\) to use the pattern command (*).

Once the per species dataset is in the laptop, I analyze them in R (See depth_species_distribution.R).

## BED with multicov information per species
Now, we need a BED file with information of coverage per species (as the sum of each sample). That's not direct, neither trivial. As argued with Enrico, the pipeline will be:
1. Get a BED_coverage per sample (run BED_coverage.sh)
2. Get BED intervals with bedops (install with conda) and use --partition. One BED interval per species
3. Join cov_column per sample and sum up
4. Callable with the sp_depth.csv (maxDepth)
5. Filter the vcf with the callable sites

### 1. Get a BED_coverage per sample
With bedtools genomecov we are calculating coverage per window in each sample, in order to laterly filter with min and max depth.

NOTE: We have to eliminate from the bamlist the samples we are not going to use, because of quality, depth or anything else:
c_ll_ca_0249
c_ll_ca_0253
c_lc_yu_0007

And also have into consideration that:
c_lr_xx_0011 is a lc

For that purpose, I moved the samples to a directory called "bad_samples" and rename c_lr_xx_0011 to c_lc_xx_0011 in the BED file (not in the BAM)

```
sample=($(ls $STORE2/lynx_genome/lynx_data/CatRef_bams/*cat_ref_sorted_rg_rmdup_sorted_indelrealigner.bam \
| rev | cut -d '/' -f 1 | rev | cut -d '_' -f 1-4))

for i in ${sample[@]}
do
echo "BED for $i"
sbatch BED_coverage.sh $i
done
```
Once we have a per sample bed, we need a per species bed (See species_bed.sh).

NOTE: In case we wanted to know how much is the genome represented in our samples bed:
```
awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}' <nombre_del_bed>
```

```
species=(lc ll lp lr)

for i in ${species[@]}
do
echo "BED for $i"
sbatch species_bed.sh $i
done
```
Although this previous command seems simple, it is presumably so computing-consuming, so we need to make it easier. For that purpose I am going to divide this work by chromosome, so we are going to get 21 files (20 chrm + rest) per species.
```
# Array of Chromosomes and species
CHR=($(cat $STORE2/reference_genomes/Felis_catus_Ref/Felis_catus.Felis_catus_9.0.dna.toplevel.fa.fai | cut -f 1 | uniq))
species=(lc ll lp lr)

# Create a BED with coordinates for each chromosome and in each species
for i in ${species[@]}
  do
  for j in ${CHR[@]:0:20}
  do
  echo "bed for $j in $i"
  sbatch chr_sp_bed2.sh $i $j
  done
  done

# For the rest of the info in the BED file:
for i in ${species[@]}
  do
  echo "bed for rest in $i"
  sbatch rest_sp_bed.sh $i
  done

```
Once we have a bed per species chromosome, we wanted to join them so that we have ONE BED PER SPECIES!
```
species=(lc ll lp lr)

#First we join each chromosome
for i in ${species[@]}
  do
  echo "bed for $i"
  sbatch sp_bed.sh $i
  done

#Then we join the rest! NOTE THAT AANG AND KZ ARE ORDERED INVERSED THAN IN REF_GENOME (Felis_catus_9)
for i in ${species[@]}
  do
  echo "Adding rest to ${i}"
  cat rest.${i}.bed >> ${i}.bed
  done
```
With this, we will have every sample ordered by chrm position but there will be overlapping regions. To control for this we are going to use bedops, which is a is a core tool for finding relationships between two or more genomic datasets. The partition operator splits all overlapping input regions into a set of disjoint segments. One or more input files may be provided; this option will segment regions from all inputs.

Because we don't have bedops in neither genomics nor CESGA server, we will use miniconda2 to create an environment where we will install it in CESGA:
```
module load miniconda2
conda create -n bedops
conda activate bedops
conda install -c bioconda bedops
```
Once created the environment called "bedops" and installed the program, we will do the partition (see bedops.sh)
With the partitioned.bed we have the intervals for every species but not the COVERAGE. Now we have to add samples coverage for each interval and then sum up the columns (total sp coverage per interval)
```
species=(lc ll lp lr)

for i in ${species[@]}
do
echo "Bedtools intersect for ${i}"
sbatch bedtools_intersect.sh $i
done        
```
Now we can assign callable vs. non callable sites. If maximum depth per species (See depth_species_distribution.R) is reached (1.5xmode), then this sites are "non-callable".
```
# define species
species=(lc ll lp lr)
# define the table
table=perspecies_depths.csv

# assign YES vs NO (callable vs non-callable) to windows based on depth
for i in ${species[@]}
 do
  max=($(grep "${i}" ${table} | cut -d',' -f5))
  echo "${i} max depth is ${max}"
  cat ${i}_partitioned.bed |
  awk -v max="${max}" '{FS="\t"; OFS="\t"; if ($4 > max) print $1, $2, $3, $4, "NO"; else print $1, $2, $3, $4, "YES";}' \
  > ${i}_partitioned_coverage_defined.bed
done

# get callable regions only (<1.5*mode), remove low-mappability and 0 depth for all samples,
# merge windows at up to 140bp distance and keep windows of at least 500 bps
for i in ${species[@]}
 do
  echo ${i}
  grep "YES" ${i}_partitioned_coverage_defined.bed |
  awk -F'\t' '{if($4 > 0) print}' |
  bedtools subtract -a - -b $STORE2/reference_genomes/Felis_catus_Ref/
  **REPETITIVE??**
  repetitive_regions/lc_rep_ALL_scaffold_coord.bed |
  bedtools merge -d 140 -i - |
  awk -F'\t' '{if($3-$2 > 500) print}' \
  > ${sp}_allsamples_LyCa_ref_chrY.callable.bed
done
```
# check amount of bp of callable regions
for i in ${species[@]}
 do
  echo ${i}
  awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}' ${sp}_allsamples_LyCa_ref_chrY.callable.bed
done

# Intersect the bed files of each species
multiIntersectBed -i $(ls *_allsamples_LyCa_ref_chrY.callable.bed) > malesamples_perspecies_LyCa_ref_chrY.callable.intersect.bed

# To check amount of bases in regions you can run:
# grep "1,2,3,4" malesamples_perspecies_LyCa_ref_chrY.callable.intersect.bed | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'
# grep "1,2,3,4" malesamples_perspecies_LyCa_ref_chrY.callable.intersect.bed | awk -F'\t' '{if($3-$2 > 500) print}' | awk -F'\t' '{print $1, $2, $3, $3-$2}' | less -S

# Generate BED file of full intersect between all species
grep "1,2,3,4" malesamples_perspecies_LyCa_ref_chrY.callable.intersect.bed > malesamples_perspecies_LyCa_ref_chrY.callable.intersect.goodinall.bed


!!!!!
 bedtools intersect -a tryA.bed -b tryB.bed -wa -wb | cut -f7 > covtryB.bed
 !!!!!
