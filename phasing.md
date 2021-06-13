
use whatshap for blocking and shapeit for haplotype phasing.
whatshap phase -o phased.vcf --reference=reference.fasta input.vcf input.bam



SHAPEIT is primarily a tool for inferring haplotypes from SNP genotypes.
It takes as input a set of genotypes and a genetic map, and produces as output, either a single set of estimated
haplotypes, or a haplotype graph that encapsulates the uncertainty about the underlying haplotypes.
