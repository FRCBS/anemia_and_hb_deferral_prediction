#!/usr/bin/env Rscript

# Read the null models for all phenotypes, and extract and represent their coefficients as a tibble.
# If no phenotypes are given on the command line, then all phenotypes will be used.

library(tidyverse)
args <- commandArgs(TRUE)

output_filename <- args[1]

#if (length(args) > 1) {
#  phenotypes <- args[2:length(args)]
#} else { 
  phenotypes <- c("hb_median", "hb_mad", "tries", "lifetime_donations", "last_two_years_donations", "S821", "P821", "P820", 
                "serious_local_ever", "fainting_ever")
#}

get_coefficients <- function(phenotype) {
  filename <- sprintf("results_%s.rda", phenotype)
  cat(sprintf("Loading from file %s\n", filename))
  tryCatch(
    error = function(cnd) { result <<- tibble(phenotype=phenotype)}, 
    {
      load(filename, verbose=FALSE)
      result <- as_tibble(t(modglmm$coefficients))
      result <- result %>% mutate(phenotype=phenotype)
    })
  return(result %>% select(phenotype, everything()))
}

df <- map_dfr(phenotypes, get_coefficients)

if (!is.na(output_filename)) {
  write_tsv(df, file = output_filename)
} else {
  print(df)  
}
