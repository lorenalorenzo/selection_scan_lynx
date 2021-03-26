## Separate by species
Create a variable with the 4 species (lc, ll, lp, lr) separated by line.
For [the variable species] in [value 0 of the variable, which is lc] do [select each sample], [select variants] **CONTINUE**
Include genotypes from this sample
```
screen -S by_species_vcf

spARRAY=($(grep -m1 "#CHR" allsamples_cat_ref.filter5.vcf | tr '\t' '\n' | grep "_" | cut -d"_" -f2 | sort -u))
for sp in ${spARRAY[@]}
  do
  samplesARRAY=($(grep -m1 "#CHROM" allsamples_cat_ref.filter5.vcf | tr '\t' '\n' | grep "${sp}"))
  /opt/gatk-4.1.0.0/gatk SelectVariants \
  -R /GRUPOS/grupolince/reference_genomes/felis_catus_genome/Felis_catus.Felis_catus_9.0.dna.toplevel.fa \
  -V allsamples_cat_ref.filter5.vcf \
  $(for j in ${samplesARRAY[@]}; do echo "-sn ${j}";done) \
  -O ${sp}_output_per_species.vcf
  done
```
### Compare if ExcessHet is recalculated
```
grep -v "#" lc_output_per_species.vcf | head -n100000 > lc_filtering_100000.vcf
grep -o -E "ExcessHet=\-?[[:digit:]]{1,8}\.?[[:digit:]]{0,8}" lc_filtering_100000.vcf  > lc_excesshet_100000.vcf
cat lc_excesshet_100000.vcf | tr '=' '\t' > lc_excesshet_toR.vcf

scp llorenzo@genomics-a.ebd.csic.es:/home/llorenzo/vcf/lc_excesshet_toR.vcf .
```
{R}
```
library(dplyr)
setwd("C:/Users/loren/Documents/R")
allsamples_excesshet <- read.table ("excesshet_toR.vcf")
lc_excesshet <- read.table ("lc_excesshet_toR.vcf")
setdiff(allsamples_excesshet, lc_excesshet)
```
NO DIFFERENCE!!
**PROBLEM: By this method, ExcessHet is not recalculated by species (nor Inbreeeding Coeff)!!**
