#!/bin/bash

# Get the common definitions
. $(dirname $0)/common.sh

# if [[ $# != 2 ]] ; then
#     input=supporting/GRCh38_positions/ALL.chr6_GRCh38.genotypes.20170504.vcf.gz
#     output=poista.koe
# else
    input=$1
    output=$2
#fi
    
$PLINK --vcf $input \
	 --vcf-half-call missing \
	 --make-bed \
	 --out $output


