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
