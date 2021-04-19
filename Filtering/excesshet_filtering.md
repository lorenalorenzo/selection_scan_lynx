---
Title: ExcessHet filtering
Author: Lorena Lorenzo Fernández
Date: 11 March, 2021
---

## ExcessHet
ExcessHet describes the heterozygosity of the called samples, giving a probability of excess heterozygosity being observed. To analyze this, we first need to separate the vcf per species (See vcf_per_species.sh).

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
With this, I obtain a vcf with this info:   ExcessHet=[value]

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
So values of ExcessHet above 13 have a p-value smaller than 0.05

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

### Compare if ExcessHet is recalculated when dividing vcf per species
```
grep -v "#" lc_output_per_species.vcf | head -n100000 > lc_filtering_100000.vcf
grep -o -E "ExcessHet=\-?[[:digit:]]{1,8}\.?[[:digit:]]{0,8}" lc_filtering_100000.vcf  > lc_excesshet_100000.vcf
cat lc_excesshet_100000.vcf | tr '=' '\t' > lc_excesshet_toR.vcf
```
Analysing in R results
```
scp llorenzo@genomics-a.ebd.csic.es:/home/llorenzo/vcf/lc_excesshet_toR.vcf .

{R}

library(dplyr)
setwd("C:/Users/loren/Documents/R")
allsamples_excesshet <- read.table ("excesshet_toR.vcf")
lc_excesshet <- read.table ("lc_excesshet_toR.vcf")
setdiff(allsamples_excesshet, lc_excesshet)
```

**PROBLEM: By this method, ExcessHet is not recalculated by species (nor Inbreeeding Coeff)!!**
