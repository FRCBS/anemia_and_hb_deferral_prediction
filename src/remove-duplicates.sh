#!/bin/bash

# Get the common definitions
. $(dirname $0)/common.sh

# removing duplicate variants

input=$1
output=$2
todrop=$3

$PLINK \
  --bfile $input \
  --exclude $todrop \
  --make-bed \
  --out $output
