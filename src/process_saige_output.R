#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(qqman))

# Input:
# - Directory of saige results
# - Phenotype name
# Output:
# - manhattan plot
# - Q-Q plot


# results_hb_median.all.summary_statistics.txt

#gwas_result_filename="saigeio/output/results_hb_median.SAIGE.vcf.genotype.txt"
args <- commandArgs(TRUE)
gwas_result_filename <- args[1]
dir <- dirname(gwas_result_filename)
phenotype <- args[2]
#gwas_result_filename <- sprintf("%s/results_%s.SAIGE.vcf.genotype.txt", dir, phenotype)
#gwas_result_filename <- sprintf("%s/results_%s.1.summary_statistics.txt", dir, phenotype)
#gwas_result_filename <- sprintf("%s/results_%s.summary_statistics.txt", dir, phenotype)
#results <- fread(gwas_result_filename, data.table=F)
results <- read_delim(gwas_result_filename, delim = " ", col_types = cols()) # Last parameter to get rid of printing of guessed col types

ymax <- ceiling(-log10(min(results$p.value))+1)

output_dir <- "../results/images"
# Manhattan plot
png(sprintf("%s/manhattan_plot_%s.png", output_dir, phenotype), height= 15, width= 35, units="cm", res=400)
manhattan(results, chr="CHR", bp="POS", p="p.value", snp="SNPID", ylim= c(0,ymax), col=c("blue", "brown"), cex = 0.6, cex.axis=0.8, las=2)
title(phenotype)
dev.off()

# http://genometoolbox.blogspot.com/2014/08/how-to-calculate-genomic-inflation.html
# https://www.mv.helsinki.fi/home/mjxpirin/GWAS_course/material/GWAS6.html
compute_gif <- function (pvalue) {
  # For p-values, calculate chi-squared statistic
  chisq <- qchisq(1 - pvalue, 1)

  # Calculate lambda gc (λgc)
  median(chisq)/qchisq(0.5, 1)
}


# Q-Q plot
png(sprintf("%s/qq_plot_%s.png", output_dir, phenotype), height= 15, width= 20, units="cm", res=300)
# Compute genomic inflation factor
lambda <- compute_gif(results$p.value)
lambda2 <- signif(lambda, 3)
qq(results$p.value, sub=bquote(~lambda == .(lambda2)))
title(phenotype)
dev.off()
