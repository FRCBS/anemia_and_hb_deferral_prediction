#!/bin/bash


set -x -e

#chromosome_files=~/FRCBS/blood_health_phewas/vcf-files.txt

#DIRFINNGEN=~/data/private/finngen
#PLINK=plink1.9

SRC=$(dirname $0)

# Get the common definitions
. $SRC/common.sh

result_dir=$DIRRESULT/pca-cleaned


input=$1
chr=$2

g=$(basename $input)
output=$result_dir/$g

#chr=$(echo $g | sed -n 's/finngen_R4_bb_\(.*\)_BLOOD.*/\1/p')
echo Chromosome is $chr

time_file=$result_dir/time-$chr.txt
> $time_file  # Make the file empty
TIME="/usr/bin/time -vao $time_file"

myprint "Cleaning data for PCA"

if false ; then 
myprint2 "Set genotype as missing if genotype probability is less than 0.95"
$TIME \
    bcftools filter -i 'FMT/GP>=0.95' --set-GTs . $input -O z -o $output
$TIME \
    bcftools index --tbi $output

myprint2 "Compute the IMPUTE2 info score"
$TIME \
    bcftools +impute-info -Ou $output | \
    bcftools query -f '%ID\t%INFO/INFO\n' > $result_dir/info-scores-$chr.txt

myprint2 "Prune variants"
$TIME \
    $PLINK \
    --vcf $output \
    --geno 0.01 \
    --maf 0.05 \
    --qual-scores $result_dir/info-scores-$chr.txt 2 1 \
    --qual-threshold 0.95 \
    --make-bed \
    --out $result_dir/pruned-chr$chr

fi

myprint2 "Finding LD-based variants"
$TIME \
    $PLINK \
    --bfile $result_dir/pruned-chr$chr \
    --indep-pairwise 500kb 50 0.1 \
    --out $result_dir/snps-for-pruning-$chr

myprint2 "Extracting LD variants"
$TIME \
    $PLINK \
    --bfile $result_dir/pruned-chr$chr \
    --extract $result_dir/snps-for-pruning-$chr.prune.in \
    --make-bed \
    --out $result_dir/ld-pruned-chr$chr

#    	--mind 0.1 \
#	--hwe 1E-6 \

#done
