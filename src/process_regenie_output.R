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

#ubuntu@blooddonorresearch02:/media/blooddonorresearch02data/R8/results$ zcat regenieio/output/step2_hb_median_hb_median.regenie.gz | head -1
#CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ N TEST BETA SE CHISQ LOG10P EXTRA

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

results <- results %>%
	mutate(P=10**-LOG10P) # LOG10P is actually -LOG10P

#ymax <- ceiling(-log10(min(results$p.value))+1) # Saige
ymax <- ceiling(max(results$LOG10P)+1) # Regenie

cat(sprintf("ymax = %f\n", ymax))

output_dir <- "../results/images"
# Manhattan plot
png(sprintf("%s/regenie_manhattan_plot_%s.png", output_dir, phenotype), height= 15, width= 35, units="cm", res=400)
manhattan(results, chr="CHROM", bp="GENPOS", p="P", snp="ID", ylim= c(0,ymax), col=c("blue", "brown"), cex = 0.6, cex.axis=0.8, las=2)
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
png(sprintf("%s/regenie_qq_plot_%s.png", output_dir, phenotype), height= 15, width= 20, units="cm", res=300)
# Compute genomic inflation factor
lambda <- compute_gif(results$P)
lambda2 <- signif(lambda, 3)
qq(results$P, sub=bquote(~lambda == .(lambda2)))
title(phenotype)
dev.off()
