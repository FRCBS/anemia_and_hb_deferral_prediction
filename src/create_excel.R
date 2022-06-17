#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(tidyverse))
#library(xlsx)
library(openxlsx)
library(lobstr)

options(warn = 1) # Show warning immediately

source("create_excel_functions.R")

args <- commandArgs(TRUE)

output <- args[1]
inputs <- args[2:length(args)]

stopifnot(str_detect(output, "\\.xlsx$"))

threshold <- 5e-8

#fast <- TRUE
fast <- FALSE

if (!fast) {
  message("Reading info scores")
  info <- get_info_scores(1:22)
  
  message("Reading enrichment data")
  enrichment <- get_enrichment()
}

dfs <- list()
for (input in inputs) {
  cat(inblue(sprintf("File %s\n", input)))
  # input <- "results.hb_median.all.summary_statistics.txt"
  phenotype <- str_match(input, "results\\.([^.]+).*")[,2]
  cat(sprintf("Phenotype %s\n", phenotype))
  df_full <- read_delim(input, delim=" ")
  print(df_full)
  n_significant <- df_full %>% filter(p.value < threshold) %>% nrow()
  n <- n_significant + 50
  df <- df_full %>% 
    arrange(p.value) %>% 
    head(n=n) %>%
    arrange(CHR, POS) %>%
    mutate(NEXT_SNPID=lead(SNPID))
  
  message("Adding annotations")
  annotations <- get_annotations(df)
  #df <- left_join(df, annotations, by=c("CHR", "POS", "REF", "ALT")) #%>%
#  df <- left_join(df, annotations, by=c("CHR", "POS", "Allele1", "Allele2")) %>%
  df <- left_join(df, annotations, by=c("CHR", "POS"))
    #select(-c(REF, ALT))

  message("Adding LDs")
  # Linkage disequilibrium between each two consecutive variants in units on R2  
  lds <- get_ld(df)
  save(lds, file=sprintf("/tmp/lds-%s.rdata", phenotype))
  df <- left_join(df, lds, by=c("SNPID","NEXT_SNPID")) %>%
    arrange(CHR, POS) %>%
    mutate(R2_TO_NEXT = ifelse(CHR != lead(CHR), 0.0, R2_TO_NEXT))  # When the next snip is in different chromosome, LD is zero
  
  if (!fast) {
    message("Adding Impute2 INFO score")
    df <- left_join(df, info, by="SNPID")
  }
  
  if (!fast) {
    df <- left_join(df, enrichment, by="SNPID")
  }

  df <- df %>% select(CHR, POS, SNPID, rsid, everything())
  
  message("Adding urls")   # Do this last so the column types of links don't get messed up
  df <- add_links(df)
  
  
  dfs[[phenotype]] <- df
  
  rm(df_full, df, annotations, lds)
  res <- gc(verbose=TRUE, full=TRUE)
  print(res)
  cat(sprintf("Memory taken by dfs: %s\n", lobstr::obj_size(dfs)))
  print(do.call(lobstr::obj_sizes, as.list(rlang::global_env())))
  #openxlsx::write.xlsx(df, file = output, overwrite=TRUE, sheetName = phenotype)
}

phenotypes <- names(dfs)
null_model_coefficients <- map2_dfr(dirname(inputs), phenotypes, get_null_model_coefficients)
dfs[["Null model coefficients"]] <- null_model_coefficients

####openxlsx::write.xlsx(dfs, file = output, overwrite=TRUE, keepNA=TRUE)
openxlsx::write.xlsx(dfs, file = output, overwrite=TRUE)
