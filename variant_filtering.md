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
ExcessHet describes the heterozygosity of the called samples, giving a probability of excess heterozygosity being observed

#See the distribution of Excess Het by analysing vcf first
```
grep -v "#" allsamples_cat_ref.filter5.vcf | less -S
```

# Number of variants with ExcessHet annotated:
FUNCIONA:
grep -v "#" allsamples_cat_ref.filter5.vcf | grep -o -E 'ExcessHet=[[:digit:]]{1,4}\.?[[:digit:]]{0,4}' | wc -l


DANI:
grep -v '#' tu_archivo.vcf | awk -F";" '{printf ("%s;%s\n", $5,$7)}' | column -t | less -S


# Extract column:
grep -v "#" allsamples_cat_ref.filter5.vcf | grep -o -E 'ExcessHet=[[:digit:]]{1,3}\.?[[:digit:]]{0,3}' | cut -d '=' -f2 > ExcessHet.table

#Distribution plot of ExcessHet in R
