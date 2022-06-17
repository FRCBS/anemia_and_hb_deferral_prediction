#!/usr/bin/env Rscript

# Find duplicate variants.

library(data.table)

args <- commandArgs(TRUE)

input <- args[1]
output <- args[2]

logfilename <- paste0(output, ".log")
# input <- './data/cleaned_genotypes/merged_cleaned_data.bim'
# output <- "./data/cleaned_genotypes/duplicates_cleaned_merged.txt"

variants <- fread(input, data.table=F)

# variants with no unique name or location:
names     <- duplicated(variants$V2) # 0
locations <- duplicated(variants[c("V1", "V4")]) # 1421 

cat(sprintf("Duplicate by name: %i, duplicate by location: %i, total: %i\n",
            sum(names), sum(locations), sum(names|locations)),
    file=logfilename)

# Write out the names of duplicates
duplicate_rows <- variants[names|locations,]
write.table(duplicate_rows$V2, file=output, sep=" ", quote=F, row.names=F, col.names=F)
