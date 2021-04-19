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
## (6) Under-represented, excessively missing variants
It's important to filter out variants which are missing completely in one or more species.
{UNDER CONSTRUCTION}

The next step should be looking for sites out of HW equilibrium (as we are interested in analysing selection).
The number of het genotypes expected under Hardy-Weinberg equilibrium is 2*(# of samples)*(ref allele frequency)*(alt allele frequency), where allele frequencies are calculated from the samples' genotypes.

## FILTERING BY HET>80% samples
First we want to extract genotype info per sample in each species.
bcftools query -f '%CHROM %POS  GTs:[ %GT]\t PLs:[ %PL]\n' lc_output_per_species.vcf > lc_genotypes
bcftools query -f '%CHROM %POS  GTs:[ %GT]\t PLs:[ %PL]\n' ll_output_per_species.vcf > ll_genotypes
