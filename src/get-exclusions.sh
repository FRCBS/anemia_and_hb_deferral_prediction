#!/bin/bash

echo -e "Chromosome\tInput variants\tInput people\tMissing people\tMissing variant\tHWE\tMAF\tOutput variants\tOutput people"


for chr in {1..22} X ; do
#for chr in 22 ; do
    input=cleaned-chr$chr.log
    input2=cleaned-deduplicated-chr$chr.log
    echo -en "$chr\t"
    sed -rn '
s/(.*) variants loaded from .bim file./\1/p
s/(.+) people \(0 males, 0 females, 20737 ambiguous\) loaded from .fam./\1/p
s/(.+) people removed due to missing genotype data \(--mind\)./\1/p
s/(.+) variants removed due to missing genotype data \(--geno\)./\1/p
s/--hwe: (.+) variants removed due to Hardy-Weinberg exact test./\1/p
s/(.+) variants removed due to minor allele threshold\(s\)/\1/p
' $input | tr "\n" "\t"
    sed -rn '
s/(.+) variants and (.*) people pass filters and QC./\1\n\2/p
' $input2 | tr "\n" "\t" | sed 's/\t$//'
    echo
done

exit




