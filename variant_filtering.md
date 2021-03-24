---
Title: Variant filtering steps for selection analysis in Lynx
Author: Lorena Lorenzo Fernández
Date: 11 March, 2021
---
Once we have our genotypes called, we want to filter the SNPs in order to get better information.
The first 5 steps of this filters are extensively explained in Enrico's Github:
https://github.com/Enricobazzi/PlanNacional_Demographic_models_Oscar/blob/master/2.Variant_Filtering.md

## (1,2,3,4,5) Applying General Filters
The script applying filters 1 through 5, which are independent of species and sequencing technology, and can therefore be applied to the whole dataset can be found at 2.Variant_Filtering-executables/General_Filter_1-2-3-4-5.sh

#### Number of variants
```
grep -v "#" allsamples_cat_ref.filter5.vcf | wc -l
22940737
````
## Separate by species
Create a variable with the 4 species (lc, ll, lp, lr) separated by line.
For [the variable species] in [value 0 of the variable, which is lc] do [select each sample], select variants **CONTINUE**
```
spARRAY=($(grep -m1 "#CHR" allsamples_cat_ref.filter5.vcf | tr '\t' '\n' | grep "_" | cut -d"_" -f2 | sort -u))
for sp in ${spARRAY[0]}
  do
  samplesARRAY=($(grep -m1 "#CHROM" allsamples_cat_ref.filter5.vcf | tr '\t' '\n' | grep "${sp}"))
  /opt/gatk-4.1.0.0/gatk SelectVariants \
  -R /GRUPOS/grupolince/reference_genomes/felis_catus_genome/Felis_catus.Felis_catus_9.0.dna.toplevel.fa \
  -V allsamples_cat_ref.filter5.vcf \
  $(for j in ${samplesARRAY[@]}; do echo "-sn ${j}";done) \
  -O ${sp}_output_per_species.vcf
  done

```
In fact, I only need to SelectVariants per species....
```
/opt/gatk-4.1.0.0/gatk IndexFeatureFile -F ll_wholegenome_LyCa_ref.sorted.filter7.vcf

/opt/gatk-4.1.0.0/gatk SelectVariants \
-R /GRUPOS/grupolince/reference_genomes/felis_catus_genome/Felis_catus.Felis_catus_9.0.dna.toplevel.fa \
-V allsamples_cat_ref.filter5.vcf \
$(for j in ${spARRAY[@]}; do echo "-sn ${j}";done) \
-O output_per_species.vcf
```

## (6) Under-represented, excessively missing variants
It's important to filter out variants which are missing completely in one or more species.
{UNDER CONSTRUCTION}

The next step should be looking for sites out of HW equilibrium (as we are interested in analysing selection).
The number of het genotypes expected under Hardy-Weinberg equilibrium is 2*(# of samples)*(ref allele frequency)*(alt allele frequency), where allele frequencies are calculated from the samples' genotypes.

## (7) Inbreeding coefficient
The output is the inbreeding coefficient 'F' (fixation) statistic, which for large sample sizes converges to the probability that an individual's two alleles are identical by descent, provided that cosanguinity is the only source of deviation from Hardy-Weinberg equilibrium.

## (8) Excess Het
ExcessHet describes the heterozygosity of the called samples, giving a probability of excess heterozygosity being observed. To analyze this, we first need to separate the vcf per species.

Previsualize vcf:
```
grep -v "#" allsamples_cat_ref.filter5.vcf | less -S
```

### Number of variants with ExcessHet annotated:
In order to know if ExcessHet is a variable present in every variant, we search for it and compare the number of lines in the original vcf with the number of lines with info about ExcessHet.
*Note that I decided that every ExcessHet value will be between 1-8 digits (ordinals and decimals) and that I typed \-? at the start due to some -0 values of ExcessHet in the vcf.
```
grep -v "#" allsamples_cat_ref.filter5.vcf | grep -o -E 'ExcessHet=\-?[[:digit:]]{1,8}\.?[[:digit:]]{0,8}' | wc -l
grep "ExcessHet" allsamples_cat_ref.filter5.vcf | wc -l
```
After this, I made a subset of 100000 variants to test the filter steps.
```
grep -v "#" allsamples_cat_ref.filter5.vcf | head -n100000 > filtering_100000.vcf
```
Then, I extract every ExcessHet value in a new vcf.
```
grep -o -E "ExcessHet=\-?[[:digit:]]{1,8}\.?[[:digit:]]{0,8}" prueba_excesshet.vcf  > excesshet_100000.vcf
```
With this, I obtain a vcf with this info:
ExcessHet=[value]

In order to get every [value] separated from the "ExcessHet" word, I have to interpret "=" as a delimiter:
```
cat excesshet_100000.vcf | tr '=' '\t' > excesshet_toR.vcf
```
This vcf is what I want to use in R to study the distribution. So, to make it easier, I am going to download this vcf in my laptop (in C:\Users\loren)
```
scp llorenzo@genomics-a.ebd.csic.es:/home/llorenzo/vcf/excesshet_toR.vcf .
```
### Distribution plot of ExcessHet in R
ExcessHet is the phred-scaled p-value for the test of heterozygosity. The phred-score is:
Phred-score= -10 * log10 (p-value)
So values of ExcessHet above 13 have a p-value smaller than 0.05.

IN GATK FILTERING I FOUND:
A] Hard-filter a large cohort callset on ExcessHet using VariantFiltration
ExcessHet filtering applies only to callsets with a large number of samples, e.g. hundreds of unrelated samples. Small cohorts should not trigger ExcessHet filtering as values should remain small. Note cohorts of consanguinous samples will inflate ExcessHet, and it is possible to limit the annotation to founders for such cohorts by providing a pedigree file during variant calling.

```
gatk --java-options "-Xmx3g -Xms3g" VariantFiltration \
    -V cohort.vcf.gz \
    --filter-expression "ExcessHet > 54.69" \
    --filter-name ExcessHet \
    -O cohort_excesshet.vcf.gz
```

This produces a VCF callset where any record with ExcessHet greater than 54.69 is filtered with the ExcessHet label in the FILTER column. The phred-scaled 54.69 corresponds to a z-score of -4.5. If a record lacks the ExcessHet annotation, it will pass filters.

Given that information, I studied the distribution of ExcessHet for the subset of 100000 variants.

{R}
```
#Upload data and packages needed for plotting
library(ggplot2)
excesshet <- read.table ("excesshet_toR.vcf")

#Plot the subset distribution
subset_dist <- ggplot(data=excesshet, aes(x=V2)) +
  geom_histogram(binwidth= 5, colour="black", fill="white") +
  scale_y_continuous(limit=c(0,70000)) +
  scale_x_continuous(breaks=c(0:170*10))

subset_dist + ggtitle("100000pb subset") + xlab("ExcessHet") + ylab("Counts")

#Plot the outliers distribution
outliers_dist <- ggplot(data=excesshet, aes(x=V2)) +
 geom_histogram(binwidth= 5, colour="black", fill="white") +
  scale_y_continuous(limit=c(0,500)) +
  scale_x_continuous(breaks=c(0:170*10))

outliers_dist + ggtitle("Outliers") + xlab("ExcessHet") + ylab("Counts")
```

Then, I decided to filter every ExcessHet value over 13 (p-value>0.05) and again look at the distribplot
```
#Plot values <13
selected_excesshet <- subset(excesshet, V2<13) #98958 of values selected (~99%)
selected_dist <- ggplot(data=selected_excesshet, aes(x=V2)) +
  geom_histogram(binwidth= 3, colour="black", fill="white") +
  scale_y_continuous(limit=c(0,70000)) +
  scale_x_continuous(breaks=c(0:170))

selected_dist + ggtitle("99%") + xlab("ExcessHet") + ylab("Counts")
```
For graphic results, go to (c:Users/loren/Documents/R)
