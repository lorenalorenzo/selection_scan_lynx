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
```
sample=($(ls $STORE2/lynx_genome/lynx_data/CatRef_bams/*cat_ref_sorted_rg_rmdup_sorted_indelrealigner.bam \
| rev | cut -d '/' -f 1 | rev | cut -d '_' -f 1-4))

for i in ${sample[@]}
do
echo "BED for $i"
sbatch BED_coverage.sh $i
done
```
Then I will have to run the previous code in an .sh format (See samtools_depth.sh). Upload the sh in the CESGA server:
```
scp /Users/lorenalorenzo/github/selection_scan_lynx/Filtering/depth/executables/BED_coverage.sh csbiellf@ft2.cesga.es:/home/csic/bie/llf
```
NOTE: We have to eliminate from the bamlist the samples we are not going to use, because of quality, depth or anything else:
c_ll_ca_0249
c_ll_ca_0253
c_lc_yu_0007

And also have into consideration that:
c_lr_xx_0011 is a lc

For that purpose, I moved the samples to a directory called "bad_samples" and rename c_lr_xx_0011 to c_lc_xx_0011

```
species=(lc ll lp lr)

cat *${species}*.coverage.bed | sort -k1,1 -k2,2n | uniq > "${species}".bed

conda create -n bedops
conda activate bedops
conda install -c bioconda bedops

for i in ${species[@]}
do
echo "BED partition for $i"
bedops --partition ${species}.bed > "${species}"_partitioned.bed
done
```
