#!/bin/bash

set -ex

SRC=$(dirname $0)

# Get the common definitions
. $SRC/common.sh

result_dir=$DIRRESULT/pca-cleaned

input_prefix=$1   # final
outliers_removed_prefix=outliers_removed_$input_prefix
pca_prefix=$outliers_removed_prefix
ibd_prefix=ibd_$outliers_removed_prefix

myprint "In file post-pca.sh"

myprint2 "Remove outliers"
$PLINK \
  --bfile $result_dir/$input_prefix \
  --remove $result_dir/outliers.tsv \
  --make-bed \
  --out $result_dir/$outliers_removed_prefix

myprint2 "Identity-by-descent"
# Generate an identity-by-descent report. It is usually best to
# perform this calculation on a marker set in approximate linkage equilibrium.
# Creates the '.genome' file
$PLINK \
  --bfile $result_dir/$outliers_removed_prefix \
  --genome \
  --min 0.2 \
  --out $result_dir/$ibd_prefix

myprint2 "Creating the 'relatives.txt' file"
$SRC/ibd_analysis.R $result_dir/$ibd_prefix.genome


myprint2 "Generate sample- and variant-based missing data reports for relatives."
# Creates the '.imiss' and '.lmiss' files.
$PLINK \
  --bfile $result_dir/$outliers_removed_prefix \
  --keep $result_dir/relatives.txt \
  --missing \
  --out $result_dir/relatives_missingness

# IBD relatives: removing individual with more missingness from each relative pair.
myprint2 "Creates the 'remove_relatives.txt' file based on ibd and missingness"
$SRC/find_bad_relatives.R $result_dir/$ibd_prefix.genome \
			  $result_dir/relatives_missingness.imiss \
			  $result_dir/$pca_prefix
# Probably the pca_prefix does not really need to be from pca result. Any bed with
# with family and individual ids will do.

$PLINK --bfile $result_dir/$outliers_removed_prefix \
		     --within $result_dir/remove_relatives.txt \
		     --pca header \
		     --pca-cluster-names nonrelative \
		     --out $result_dir/ultimate-pca
