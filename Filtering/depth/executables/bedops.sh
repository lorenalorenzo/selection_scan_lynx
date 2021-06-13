#!/bin/bash
#SBATCH -t 10:00:00
#SBATCH -c 1

#We want to create unsolapped intervals in each per species bed, in other words, join intervals with the same coverage.

#First we need to load the environment where we have bedops installed and activate it.
module load cesga/2018
module load miniconda2
conda activate bedops

#Calculate intervals with partition option of bedops
species=(lc ll lp lr)

for i in ${species[@]}
do
echo "BED partition for $i"
bedops --partition $LUSTRE/${i}.bed > $LUSTRE/"${i}"_partitioned.bed
done
