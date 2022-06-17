#!/bin/bash

# Convert the output of Saige to a form that can be fed to VEP of Ensembl.

# Format (basically VCF):
# 1 65568 . A C . . .
# 2 265023 . C T . . .
# 3 319780 . GA G . . .

input=$1
echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
awk '{printf "%i\t%i\t.\t%s\t%s\t.\t.\t.\n", $1, $2, $4, $5}' $input
