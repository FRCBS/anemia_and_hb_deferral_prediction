#!/usr/bin/env Rscript

#library(GADMTools)
#library(gridExtra)
#library(ggpubr)
#library(tidyverse)
#library(viridis)
library(knitr)
#library(RColorBrewer)
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(ggpubr))

args <- commandArgs(TRUE)
prefix <- args[1]
dir_1kg <- args[2]
outliers_filename <- args[3]

#prefix <- "/home/toivoja/FRCBS/blood_health_phewas/results/pca-cleaned/pca-FG-1KG-intersection"
#dir_1kg <- "~/data/public/1000genomes"

# This contains the PCA result
result_dir = dirname(prefix)
eigenvector_filename <- paste(prefix, ".eigenvec", sep="")
eigenvalue_filename  <- paste(prefix, ".eigenval", sep="")
aberrant_outliers_filename    <- paste(prefix, "aberrant_pca_outliers.png", sep="_")
aberrant_results_filename     <- paste(prefix, "aberrant_results.rdata", sep="_")
#outliers_filename <- paste(result_dir, "outliers.tsv", sep="/")

recompute_results <- TRUE

pca <- read_delim(eigenvector_filename, delim=" ", col_names = TRUE)

# Get population labels from 1000 Genomes data
population_filename <- paste(dir_1kg, "integrated_call_samples_v3.20130502.ALL.panel", sep="/")
# The input file has stupidly two tab characters at the end of the header
population <- read_delim(population_filename, delim="\t", col_names = c("sample", "pop", "super_pop", "gender"))


# Create population name for each individual.
# Use super population except for Finns use the population (FIN).
population <- population %>% 
  mutate(super_pop=ifelse(pop=="FIN", "FIN", super_pop)) %>%
  select("IID"="sample", "population"="super_pop")

# Create population labels for individuals (samples) in the combined FinnGen and 1KG data sets
ids <- tibble(IID=pca$IID)
population <- left_join(ids, population, by="IID") %>%
  replace_na(list(population="FinnGen")) # Those not in 1kg are in FinnGen

# Join the IID/population dataframe and the PCA dataframes
pca_and_population <- pca %>% inner_join(population, by="IID")


population_names <- unique(c("FIN", "FinnGen", pca_and_population$population))
#temp <- tibble(population=population_names)
#pca <- inner_join(temp, pca, by="population")  # Just to reorder
pca_and_population <- pca_and_population %>% mutate(population = factor(pca_and_population$population, levels=population_names))
#pca <- pca %>% select(everything(), "population")

# This order defines the plotting order
pca_and_population <- arrange(pca_and_population, desc(population))
#print(head(pca))
#summary(pca)

# Pairs of the first five principal components to be plotted
#column.combinations <- combn(1:5, 2)  # Creates a 2x10 matrix of numbers 1,...,5
column.combinations <- combn(1:2, 2)

# Plot the ith pair of principal components
aberrant_function <- function(i){
  xi <- column.combinations[1,i]
  yi <- column.combinations[2,i]
  xcolname <- sprintf("PC%i", xi)
  ycolname <- sprintf("PC%i", yi)
  x <- pca_and_population %>% filter(population %in% c("FinnGen")) %>% dplyr::select(all_of(xcolname), all_of(ycolname))
  res <- aberrant::aberrant (x, lambda=30, alpha=1, beta=20, standardize= FALSE, uncorr = TRUE )
  res
}

if (recompute_results || !file.exists(aberrant_results_filename)) {
  #future::plan(future::multisession, workers = 4)
  #results <- furrr::future_map(1:ncol(column.combinations), aberrant_function, .options = furrr::furrr_options(seed = TRUE))
  results <- map(1:ncol(column.combinations), aberrant_function)
  save(results, file=aberrant_results_filename)
} else {
  load(aberrant_results_filename, verbose=TRUE)
}

# Visualize the results. The ellipse is the 99% confidence area.

# Making a combined plot does not seem to work as aberrant.plot does not return anything.

#aberrant_plots <- map(results, aberrant::aberrant.plot)
#combined_aberrant_plot <- ggpubr::ggarrange(plotlist=aberrant_plots, common.legend=T, ncol=2, nrow=5, legend="bottom", labels = paste(LETTERS[1:length(aberrant_plots)], ")", sep=""))

png(aberrant_outliers_filename, height= 300, width= 300, units="px")
aberrant::aberrant.plot(results[[1]])
dev.off()

#ggsave(aberrant_outliers_filename, device="png", height= 40, width= 30, units="cm")



# We use the first two principal components to decide the outliers.
# Verify that we get the same number of outliers as in the figure.


group <- results[[1]]$group
table(group)

# Get the ids of the outliers.

outliers <- pca_and_population %>% filter(population %in% c("FinnGen")) %>% filter(group==1)
write_tsv(outliers %>% select(FID, IID), outliers_filename, col_names=FALSE)
cat(sprintf("Wrote outliers to file %s\n", outliers_filename))


