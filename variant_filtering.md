---
Title: Variant filtering steps for selection analysis in Lynx
Author: Lorena Lorenzo Fernández
Date: 11 March, 2021
---
# VARIANT FILTERING
Once we have our genotypes called, we want to filter the SNPs in order to get better information.
The first 5 steps of this filters are extensively explain in Enrico's Github:
https://github.com/Enricobazzi/PlanNacional_Demographic_models_Oscar/blob/master/2.Variant_Filtering.md

## Workflow
### (1,2,3,4,5) Applying General Filters
The script applying filters 1 through 5, which are independent of species and sequencing technology, and can therefore be applied to the whole dataset can be found at 2.Variant_Filtering-executables/General_Filter_1-2-3-4-5.sh

#### Number of variants
```
grep -v "#" allsamples_cat_ref.filter5.vcf | wc -l
22940737
````

### (6) Under-represented, excessively missing variants
It's important to filter out variants which are missing completely in one or more species.
{UNDER CONSTRUCTION}

The next step should be looking for sites out of HW equilibrium (as we are interested in analysing selection).
The number of het genotypes expected under Hardy-Weinberg equilibrium is 2*(# of samples)*(ref allele frequency)*(alt allele frequency), where allele frequencies are calculated from the samples' genotypes.

### (7) Inbreeding coefficient
The output is the inbreeding coefficient 'F' (fixation) statistic, which for large sample sizes converges to the probability that an individual's two alleles are identical by descent, provided that cosanguinity is the only source of deviation from Hardy-Weinberg equilibrium.

### (8) Excess Het
ExcessHet describes the heterozygosity of the called samples, giving a probability of excess heterozygosity being observed. To analyze this, we first need to separate the vcf per species.


#See the distribution of Excess Het by analysing vcf first
```
grep -v "#" allsamples_cat_ref.filter5.vcf | less -S
```

# Number of variants with ExcessHet annotated:
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
# Distribution plot of ExcessHet in R
