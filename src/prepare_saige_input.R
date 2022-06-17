#!/usr/bin/env Rscript

library(data.table)
suppressPackageStartupMessages(library(tidyverse))
library(parallel)
library(gridExtra)
library(ggpubr)
library(qqman)
library(SPAtest)

# Input:
# - eigenvec file (for principal components)
# - phenotypes.txt
# - covariates.txt (non-principal components)
# Output:
# - the combined phenotype-covariate file

args <- commandArgs(TRUE)

#result_dir <- '~/FRCBS/blood_health_phewas/results'


pca_filename       <- args[1] # final_pca_ld_pruned_cleaned_merged_1_23.eigenvec
phenotype_filename <- args[2]
covariate_filename <- args[3]
output_filename    <- args[4]
#phenotype_filename <- paste(result_dir, "phenotypes.txt", sep="/")
#covariate_filename <- paste(result_dir, "covariates.txt", sep="/")
#output_filename <- paste(result_dir, "saigeio/input/pheno_covar.tsv", sep="/")

pca <- read_delim(pca_filename, delim=" ", col_names = TRUE)
print(pca)
pca <- pca %>% mutate(donor=IID)



phenotypes <- read_delim(phenotype_filename, delim="\t", col_names = TRUE)
covariates <- read_delim(covariate_filename, delim="\t", col_names = TRUE)
print(phenotypes)
print(covariates)
phenotypes <- phenotypes %>% drop_na(hb_median)   # Phenotypes cannot be NA

phenotypes_donors <- phenotypes$donor
pca_donors <- pca$donor
donors <- intersect(pca_donors, phenotypes_donors)

phenotypes <- phenotypes %>% 
  filter(donor %in% donors) 
pca <- pca %>% 
  filter(donor %in% donors)
covariates <- covariates %>% 
  filter(donor %in% donors)

pca_and_covariates <- left_join(pca, covariates, by="donor") %>%
  select(c("donor", sprintf("PC%i", 1:10), "sex", "age", "weight", "height", "smoking"))
# Covariates cannot have NA values
does_column_contain_NA <- map_lgl(pca_and_covariates, function(v) any(is.na(v)))
print(does_column_contain_NA)
stopifnot(any(does_column_contain_NA) == FALSE)

#phenotypes <- phenotypes %>% arrange(donor)
#pca_and_covariates <- pca_and_covariates %>% arrange(donor)

pheno_covar <- inner_join(phenotypes, pca_and_covariates) %>%
  rename(IID=donor)
#print(pheno_covar)

write_delim(pheno_covar, output_filename, delim="\t", col_names = TRUE)




