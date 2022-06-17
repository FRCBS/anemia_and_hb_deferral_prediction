#!/usr/bin/env Rscript

library(data.table)
library(tidyverse)
library(parallel)

# Input:
# - ibd_ld_pruned_cleaned_merged_1_23.genome
# Output:
# - relatives.txt

#input <- '/media/tk-tutk/TOSHIBA_3T_BU01/verenluovuttajat/pca/ibd_ld_pruned_cleaned_merged_1_23.genome'
args <- commandArgs(TRUE)
input_filename <- args[1]
result_dir <- dirname(input_filename)
output_filename <- paste(result_dir, "relatives.txt", sep="/")

#############################################################################################################################################
# IBD analysis
relatives <- fread(input_filename, data.table=F)
# first degree relatives for removal: PI_HAT value cut-off 0.4
first_degree <- relatives[relatives$PI_HAT >= 0.375,]

first_subject <- first_degree[,1:2]
second_subject <- first_degree[,3:4]
colnames(second_subject) <- colnames(first_subject)

all_related <- rbind(first_subject, second_subject)
write.table(all_related, file=output_filename, sep="\t", quote=F, row.names=F, col.names=F)



