#!/bin/bash

set -ex

covars=PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,sex,age,weight,height,smoking

#phenotypes="hb_median hb_mad tries lifetime_donations last_two_years_donations S821 P821 serious_local_ever fainting_ever better_fainting_ever bin_S821 bin_P821 bin_P820" # P820

#phenotypes="better_fainting_ever"
phenotypes="$@"

binary_phenotypes="serious_local_ever fainting_ever better_fainting_ever bin_S821 bin_P821 bin_P820 bin_S821_or_P821"



#for pheno in ${phenotypes//,/ } ; do 
for pheno in $phenotypes ; do 
#for pheno in hb_median ; do 

    echo Phenotype is $pheno
    
    output_prefix=./output/results.$pheno

    if grep $pheno <<< $binary_phenotypes ; then
    	echo Binary phenotype
    	#For Binary traits:
    	Rscript step1_fitNULLGLMM.R     \
    		--plinkFile=./input/input \
    		--phenoFile=./input/pheno_covar.tsv \
    		--phenoCol=$pheno \
    		--covarColList=$covars \
    		--sampleIDColinphenoFile=IID \
    		--traitType=binary        \
    		--outputPrefix=$output_prefix \
    		--nThreads=20 \
    		--LOCO=FALSE \
    		--traceCVcutoff 0.0025 \
    		--ratioCVcutoff 0.001 \
    		--IsOverwriteVarianceRatioFile=TRUE ## v0.38. Whether to overwrite the variance ratio file if the file already exists

    else
    	echo Quantitative phenotype
    	# For Quantitative traits, if not normally distributed, inverse normalization
    	# needs to be specified to be TRUE --invNormalize=TRUE
    	Rscript step1_fitNULLGLMM.R     \
    		--plinkFile=./input/input \
    		--phenoFile=./input/pheno_covar.tsv \
    		--phenoCol=$pheno \
    		--covarColList=$covars \
    		--sampleIDColinphenoFile=IID \
    		--traitType=quantitative       \
    		--invNormalize=TRUE	\
    		--outputPrefix=$output_prefix\
    		--nThreads=20 \
    		--LOCO=FALSE	\
    		--traceCVcutoff 0.0025 \
    		--ratioCVcutoff 0.001 \
    		--tauInit=1,0 \
    		--IsOverwriteVarianceRatioFile=TRUE ## v0.38. Whether to overwrite the variance ratio file if the file already exists
    fi

    echo -e "\n\nStep 2\n\n"

    chromosomes=$(seq 1 23)

    # Note: the X chromosome is 23
    for chromosome in $chromosomes ; do
    #VCF containing genotypes (--vcfField=GT)
	Rscript step2_SPAtests.R \
            --vcfFile=./input/genotypes-chr$chromosome.vcf.gz \
            --vcfFileIndex=./input/genotypes-chr$chromosome.vcf.gz.tbi \
            --vcfField=GT \
	    --chrom=$chromosome \
	    --minMAF=0.000001 \
            --minMAC=5 \
            --GMMATmodelFile=$output_prefix.rda \
            --varianceRatioFile=$output_prefix.varianceRatio.txt \
            --SAIGEOutputFile=${output_prefix}.$chromosome.summary_statistics.txt \
            --numLinesOutput=1000 \
            --IsOutputAFinCaseCtrl=TRUE \
	    --IsOutputlogPforSingle=FALSE \
	    --analysisType=additive

    #        --chrom=1 \
	#        --minMAF=0.0001 \
	#        --minMAC=1 \
	#        --sampleFile=./input/sampleIDindosage.txt \

    done
done
