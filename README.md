# selection_scan_lynx
In this chapter I am, by different ways and steps, trying to analyze positive selection in lynx.

Here I am going to copy some information that could be useful for myself in a while

*Daniel Kleinman told me to use awk which seems to be impossible if each variant have diferent information (thus, column 6 may indicate eg. allele frequency in one variant while could be inbreeding coefficient in another.)
grep -v '#' tu_archivo.vcf | awk -F";" '{printf ("%s;%s\n", $5,$7)}' | column -t | less -S
grep -v "#" allsamples_cat_ref.filter5.vcf | grep -o -E 'ExcessHet=[[:digit:]]{1,3}\.?[[:digit:]]{0,3}' | cut -d '=' -f2 > ExcessHet.table
