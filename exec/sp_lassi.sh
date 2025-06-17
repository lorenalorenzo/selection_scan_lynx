#!/bin/bash
#SBATCH --time=10:00:00
#SBATCH --mem=10GB

#name variables
sp=($(echo $1))

lassip \
       --spectra $LUSTRE/selection_scan/saltiLASSI/*${sp}.lassip.hap.spectra.gz \
       --salti \
       --out $LUSTRE/selection_scan/saltiLASSI/${sp}_salti
