##########################################
######## Separate a VCF by species #######
##########################################

# Create a variable with the 4 species (lc, ll, lp, lr) separated by line.
# For [the variable species] in [value 0 of the variable, which is lc] do [select each sample], [select variants]
# Include genotypes from this sample

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

REVISAR ESTO, LO QUE YO HE HECHO ES LO DE ARRIBA, AHORA TENGO QUE ADAPTARLO

  for i in lc ll lp lr
   do
    echo ${i}
    samples=($(ls $LUSTRE/samples_subset_depth/*masked.depth |rev | cut -d '/' -f 1 | grep "${i}" | rev | cut -d '.' -f 1)))
    /opt/gatk-4.1.0.0/gatk SelectVariants \
    -R /GRUPOS/grupolince/reference_genomes/felis_catus_genome/Felis_catus.Felis_catus_9.0.dna.toplevel.fa \
    -V allsamples_cat_ref.filter5.vcf \
    $(for j in ${samplesARRAY[@]}; do echo "-sn ${j}";done) \
    -O ${sp}_output_per_species.vcf
    done
