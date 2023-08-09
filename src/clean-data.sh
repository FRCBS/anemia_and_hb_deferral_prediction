#!/bin/bash

SRC=$(dirname $0)

# Get the common definitions
. $SRC/common.sh

input=$1
output=$2
output2=$3
chromosome=$4

#DIRRESULT=$(dirname $output)

if [[ $chromosome == "X" ]] ; then
    HWE=""
else
    HWE="--hwe 1E-6"
fi

$PLINK \
    --bfile $input \
    --geno 0.05 \
    --mind 0.1 \
    $HWE \
    --maf 0.000001 \
    --mac 5 \
    --make-bed \
    --out $output

#    --maf 0.01 \


# Remove multialleles
input2=$output
#output2=$DIRRESULT/cleaned-deduplicated-chr$chromosome
$SRC/find_duplicates.R $input2.bim $DIRRESULT/duplicates_to_drop-chr$chromosome.txt
$SRC/remove-duplicates.sh $input2 $output2 $DIRRESULT/duplicates_to_drop-chr$chromosome.txt

rm ${output}.bed
