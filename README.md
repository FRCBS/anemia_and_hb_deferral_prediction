Source code of the manuscript "The value of genetic data from 665 460 individuals in predicting anemia and suitability to donate blood".

# Description of files

## Top-level

### ./

File | Description
-----|------------
pipeline-phewas-30k.svg | Description of the R5 pipeline

## Exploratory analysis prior to running GWAS with Saige

### src/

File | Description
-----|------------
anemia_metaresults_analysis.Rmd Explore the results of the FG+UKBB meta-analysis of IDA
compare_progesa_and_finngen.Rmd
compare_finngen_and_1kg.Rmd
explore_deferral_dataframe.Rmd | Explore the contents of the deferral dataframe
explore_haitat_dataframes.Rmd | Explore the contents of the defects dataframe
venn_diagram.Rmd | Simple test of using the RVenn package
compute_phenotypes.Rmd | Compute the phenotypes for the GWAS

## Saige pipeline

### src/

File | Description
-----|------------
gwas_pipeline.sh | Runs the full Saige pipeline
common.R | Common functions an dataframe used in many places
common.sh | Contains the paths of the 1KG, FinnGen and result files
epouta-common.sh | Contains the paths on the ePouta environment
do-pca-and-grm-data.sh | Pipeline for PCA and GRM
clean-data-for-pca.sh | Perform quality control for PCA and GRM
find_duplicates.R | Create a set of duplicate variants
find_outliers.R | Find and visualize the PCA outliers
post-pca.sh | Pipeline that removes outliers, does IBD analysis and performs PCA for the second time
ibd_analysis.R | Perform the Independence by Descent analysis
find_bad_relatives.R | Create a set of individuals to be removed, so that no too close relatives exists
convert-to-bed.sh | Convert VCF files to BED format
clean-data.sh | Perform quality control for GWAS
remove-duplicates.sh | Remove duplicate variants
prepare_saige_input.R | Combine that PCA, covariate and phenotype dataframes
process_saige_output.R | Create the Manhattan and QQ plots after Saige has run

### src/saigeio/

File | Description
-----|------------
run.sh | Runs both steps of Saige
step1_fitNULLGLMM.R | Saige step1
step2_SPAtests.R | Saige step2

## Analysis done after Saige pipeline

### src/

File | Description
-----|------------
anemia_logistic_and_cox.Rmd | Fit logistic regression and Cox models for iron deficiency anemia
deferral_logistic_and_cox.Rmd | Fit logistic regression and Cox models for Hb deferral
pca-on-map-of-finland.Rmd | Plot the PCA results and outliers.
FIN_adm4.sf.rds | Data of Finnish geographical map
prepare-for-vep.sh | Convert Saige summary stats to a form that VEP requires
get_coefficients.R | Extract the coefficients of the Saige null models
get-exclusions2.sh | Extract the data exclusion counts from the PLINK log files
get-exclusions.sh | Another version of the above
create_excel_functions | R Helper functions
create_excel.R | Create an annotated Excel files of the GWAS results
create_anemia_and_deferral_article_results.Rmd | Create tables and figures for the anemia article

