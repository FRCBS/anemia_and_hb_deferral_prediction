#!/bin/bash

# Stop immediately on error (-e). Echo each command line (-x).
set -ex

SRC=$(dirname $0)

# Get the common definitions
. $SRC/common.sh

result_dir=$DIRRESULT/pca-cleaned

myprint "In file do-pca-and-grm-data.sh"

# Clean and prune FinnGen data
parallel -j20 $SRC/clean-data-for-pca.sh \
	 $DIRFINNGEN/BB_finngen_R7_chr{}_BLOOD_SERVICE_BB_extracted_consentdonors.vcf.gz \
	 {} \
	 ::: {1..22}
#	 $DIRFINNGEN/finngen_R4_bb_chr{}_BLOOD_SERVICE_extracted_*_donor.vcf.gz \



myprint2 "Concatenate FinnGen chromosomes"
for chr in {1..22} ; do
    echo $result_dir/ld-pruned-chr$chr
done > $result_dir/files-to-join.txt
$PLINK \
    --make-bed \
    --merge-list $result_dir/files-to-join.txt \
    --out $result_dir/final

myprint2 "Convert to vcf since merging with 1KG did not work using plink and bed files"
$PLINK --bfile $result_dir/final \
       --recode vcf-iid bgz \
       --out $result_dir/final-vcf
bcftools index --tbi $result_dir/final-vcf.vcf.gz

myprint2 "Concatenate 1KG chromosomes"
bcftools concat -O z -o $DIR1KG/1KG-first-22-chromosomes.vcf.gz \
	 $DIR1KG/ALL.chr{1..22}_GRCh38.genotypes.20170504.vcf.gz
bcftools index --tbi $DIR1KG/1KG-first-22-chromosomes.vcf.gz

myprint2 "Find common variants. Identifying keys are: CHROM, POS, REF, and ALT:"
bcftools isec -c none -n=2  \
	 $result_dir/final-vcf.vcf.gz \
	 $DIR1KG/1KG-first-22-chromosomes.vcf.gz > $result_dir/common-positions.txt

myprint2 "Merge FinnGen and 1KG to common variants"
bcftools merge -m all -R $result_dir/common-positions.txt \
	 -O z -o $result_dir/FG-1KG-intersection.vcf.gz \
	 $result_dir/final-vcf.vcf.gz \
	 $DIR1KG/1KG-first-22-chromosomes.vcf.gz
bcftools index --tbi $result_dir/FG-1KG-intersection.vcf.gz

myprint2 "Do the PCA for the combination of FG and 1KG"
$PLINK							\
    --vcf $result_dir/FG-1KG-intersection.vcf.gz	\
    --pca header					\
    --out $result_dir/pca-FG-1KG-intersection



myprint2 "Find outliers"
$SRC/find_outliers.R $result_dir/pca-FG-1KG-intersection \
		     $DIR1KG/../../ \
		     $result_dir/outliers.tsv

myprint2 "Remove outliers and redo PCA"
$SRC/post-pca.sh final
