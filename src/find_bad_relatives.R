#!/usr/bin/env Rscript

library(data.table)
library(tidyverse)
library(parallel)

# Input:
# - ibd_ld_pruned_cleaned_merged_1_23.genome
# - relatives_missingness.imiss
# Output:
# - remove_relatives.txt

args <- commandArgs(TRUE)

# input <- 'ibd_ld_pruned_cleaned_merged_1_23.genome'
# missingness_filename <- 'relatives_missingness.imiss'
input <- args[1]  # IBD
missingness_filename <- args[2]
pca_prefix <- args[3]
result_dir <- dirname(input)
output_filename <- paste(result_dir, "remove_relatives.txt", sep="/")

relatives <- fread(input, data.table=F)
first_degree <- relatives[relatives$PI_HAT >= 0.375,] # First-degree relatives are either siblings or in child-parent relationship
# missingness for individuals with ./src/pca_and_ibd.sh
missingness <- fread(missingness_filename, data.table=F)

# Between first-degree relatives, drop those with higher missing genotype rate
individuals_for_removal <- function(x){
  individual_1 <- missingness[missingness$IID == first_degree[x,2],]
  individual_2 <- missingness[missingness$IID == first_degree[x,4],]
  if(individual_2$F_MISS > individual_1$F_MISS){
    return(individual_2$IID)
  } else {
    return(individual_1$IID)
  }
}
# Iterate over all first-degree relative pairs
remove_these <- sapply(1:nrow(first_degree), individuals_for_removal)  
remove_these <- unique(remove_these)

# The PCA input is used here only as a reference: either to restrict the set of individuals
# or to get both FID and IID.
ids_filename <- paste(pca_prefix, ".fam", sep="")
#pca <- fread(eigenvector_filename, data.table=F)
fam <- read_delim(ids_filename, col_names=FALSE, delim=" ")
colnames(fam) <- c("FID", "IID", "X3", "X4", "X5", "X6")
#rows <- fam$IID %in% remove_these
#list_for_removal <- pca[rows,1:2]
df_for_removal <- fam %>%
	       mutate(cluster = ifelse(IID %in% remove_these, "relative", "nonrelative")) %>%
	       select(FID, IID, cluster)
write_delim(df_for_removal, output_filename, delim="\t", col_names=F)




