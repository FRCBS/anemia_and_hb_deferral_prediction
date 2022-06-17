#!/bin/bash

# Stop immediately on error (-e). Echo each command line (-x).
set -ex

SRC=$(dirname $0)

# Get the common definitions
. $SRC/common.sh

mkdir -p $DIRRESULT/{pca-cleaned,images,saigeio/input,saigeio/output}

#RERUN_PIPELINE=1
unset RERUN_PIPELINE    # Just run Saige

if [ -v RERUN_PIPELINE ] ; then   # If set rerun everything

# Preprocess data for GRM and PCA, and perform PCA
$SRC/do-pca-and-grm-data.sh



myprint2 "Convert FinnGen VCFs to BEDs"
for chromosome in {1..23} ; do
    input=$DIRFINNGEN/BB_finngen_R7_chr${chromosome}_BLOOD_SERVICE_BB_extracted_consentdonors.vcf.gz
    #input=$(ls $DIRFINNGEN/finngen_R4_bb_chr${chromosome}_BLOOD_SERVICE_extracted_*_donor.vcf.gz)
    output=$DIRRESULT/initial-chr$chromosome
    sem -j 20 $SRC/convert-to-bed.sh $input $output
done

sem --wait




myprint2 "Clean Finngen data"
for chromosome in {1..23} ; do
    input=$DIRRESULT/initial-chr$chromosome
    output=$DIRRESULT/cleaned-chr$chromosome
    output2=$DIRRESULT/cleaned-deduplicated-chr$chromosome
    sem -j 20 $SRC/clean-data.sh $input $output $output2 $chromosome
done

sem --wait


fi


myprint "Prepare input for SAIGE"

myprint2 "Create pheno_covar file."
# Its name is the last filename on the below command line
$SRC/prepare_saige_input.R $DIRRESULT/pca-cleaned/ultimate-pca.eigenvec \
			   $DIRRESULT/phenotypes.txt \
			   $DIRRESULT/covariates.txt \
			   $DIRRESULT/saigeio/input/pheno_covar.tsv

if [ -v RERUN_PIPELINE ] ; then   # If set rerun everything

myprint2 "Copy input for step 1 of Saige"
for end in bed bim fam ; do 
    cp $DIRRESULT/pca-cleaned/outliers_removed_final.$end \
       $DIRRESULT/saigeio/input/input.$end
done



myprint2 "Copy input for step 2 of Saige"
# Create input in vcf form
parallel -j 20 $PLINK --bfile $DIRRESULT/cleaned-deduplicated-chr{} \
	 --recode vcf-iid bgz \
	 --out $DIRRESULT/saigeio/input/genotypes-chr{} :::  {1..23}
#mv $DIRRESULT/saigeio/input/genotypes-chrX.vcf.gz $DIRRESULT/saigeio/input/genotypes-chr23.vcf.gz

parallel -j 20 bcftools index --tbi $DIRRESULT/saigeio/input/genotypes-chr{}.vcf.gz :::  {1..23}

fi


myprint2 "copy the needed scripts from git repository to saige directory"
cp $SRC/saigeio/*.R $SRC/saigeio/run.sh $DIRRESULT/saigeio

#phenotypes="hb_median hb_mad tries lifetime_donations last_two_years_donations S821 P821 serious_local_ever fainting_ever better_fainting_ever bin_S821 bin_P821 bin_P820 bin_S821_or_P821" # P820
#phenotypes="S821"
#phenotypes="bin_S821 bin_P821 bin_P820"
#phenotypes="better_fainting_ever"
phenotypes="bin_S821_or_P821"

myprint2 "Run the GWASes."
#Share the $DIRRESULT/saigeio directory between the host machine and the container
docker container run -v $DIRRESULT/saigeio:/saigeio --rm -w /saigeio finngen/saige:0.39.1.fg bash run.sh $phenotypes


# Combine results over chromosomes
for phenotype in $phenotypes ; do 
	sed '1p;/^CHR/d' $DIRRESULT/saigeio/output/results.$phenotype.{1..23}.summary_statistics.txt \
	> $DIRRESULT/saigeio/output/results.$phenotype.all.summary_statistics.txt
done

myprint2 "Draw the Manhattan and QQ plots for each phenotype"

parallel -j 7 $SRC/process_saige_output.R \
	 $DIRRESULT/saigeio/output/results.{}.all.summary_statistics.txt {} ::: $phenotypes

