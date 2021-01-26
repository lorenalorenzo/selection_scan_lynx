---
Title: First steps analysing selection in Lynx
Author: Lorena Lorenzo Fernández
Date: 25 January, 2021
---
### Creating a folder in the terminal and seeing the content
```
mkdir vcf
ls
```

### Copying the data from Enrico home to my pc and upload to my folder
```
scp ebazzicalupo@genomics-a.ebd.csic.es:/home/ebazzicalupo/LC_LR/CatRef_vcfs/lc_cat_ref.filter5.vcf .
lc_cat_ref.filter5.vcf                                                  
grep "#" lc_cat_ref.filter5.vcf > prueba
grep -v "#" lc_cat_ref.filter5.vcf | head -10000 >> prueba
less -S prueba
grep -v "AC=0;" > prueba_noac0
^C
grep -v "AC=0;" prueba > prueba_noac0
wc -l prueba
14552 prueba
wc -l prueba_noac0
6852 prueba_noac0
grep -v "AF=1;" prueba_noac0 | less -S
grep "AF=1;" prueba_noac0 | less -S
grep "AF=1;" prueba_noac0 | less -S
less -S prueba_noac0
mv prueba_noac0 prueba_noac0.vcf
```
