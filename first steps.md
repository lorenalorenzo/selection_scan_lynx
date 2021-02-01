---
Title: First steps analysing selection in Lynx
Author: Lorena Lorenzo Fernández
Date: 25 January, 2021
---
## PREPARING THE DATA
### Creating a folder in the terminal and seeing the content
```
mkdir vcf
ls
```

### Copying the data from Enrico home to my pc and upload to my folder
```
scp ebazzicalupo@genomics-a.ebd.csic.es:/home/ebazzicalupo/LC_LR/CatRef_vcfs/lc_cat_ref.filter5.vcf .
lc_cat_ref.filter5.vcf     
```
### Subsetting the data for trials                                             
```
grep "#" lc_cat_ref.filter5.vcf > prueba
grep -v "#" lc_cat_ref.filter5.vcf | head -10000 >> prueba
grep -v "AC=0;" prueba > prueba_noac0
```
This results in:
wc -l prueba
14552 prueba
wc -l prueba_noac0
6852 prueba_noac0

### Change the name to make a vcf file
```
mv prueba_noac0 prueba_noac0.vcf
```
## PREPARING THE SCRIPT (rehh)
```
install.packages("vcfR")
install.packages("rehh")
library(vcfR)
library(rehh)
trial1<- data2haplohh(hap_file= "prueba_noac0.vcf", polarize_vcf= FALSE, vcf_reader= "vcfR")
```
OUTPUT
* Reading input file(s) *
Using package 'vcfR' to read vcf.
Extracting map information.
Scanning file to determine attributes.
File attributes:
  meta lines: 4551
  header_line: 4552
  variant count: 2300
  column count: 29
Meta line 4551 read in.
All meta lines processed.
gt matrix initialized.
Character matrix gt created.
  Character matrix gt rows: 2300
  Character matrix gt cols: 29
  skip: 0
  nrows: 2300
  row_num: 0
Processed variant: 2300
All variants processed
Extracting haplotypes.
Number of individuals which are
Haploid Diploid Triploid, ... :
1 2
0 20
No marker identifiers found in vcf file.
* Filtering data *
Discard markers genotyped on less than 100 % of haplotypes.
517 markers discarded.
1783 markers remaining.
Data consists of 40 haplotypes and 1783 markers.
Number of mono-, bi-, multi-allelic markers:
1 2
0 1783

```
scan<- scan_hh(trial1, polarized= FALSE)
ihs<- ihh2ihs(scan, freqbin= 1)

plot<- freqbinplot(ihs)
plot_ihs<-distribplot(ihs$ihs$IHS, xlab="iHS")

dev.print(postscript,file="distribution.ihs")
