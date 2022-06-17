inblue <- function(s) {
  return(paste0("\U1b[36m", s, "\U1b[0m"))
}

get_info_scores <- function(chromosomes) {
  dir <- "~/FRCBS/blood_health_phewas/results/pca-cleaned" 
  df <- map_dfr(chromosomes, function(i) read_tsv(sprintf("%s/info-scores-chr%i.txt", dir, i), col_names = FALSE))
  colnames(df) <- c("SNPID", "INFO")
  return(df)
}

# Get annotations using Variant Effect Predictor
# The input parameter 'df' should have columns CHR, POS, ID, Allele1, Allele2
get_annotations <- function(df) {
  #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO
  input <- df %>% 
    mutate(ID=".", QUAL=".", FILTER=".", INFO=".") %>%
    select("#CHROM"="CHR", POS, ID, REF=Allele1, ALT=Allele2, QUAL, FILTER, INFO)
  #input_filename <- "/tmp/input.vcf"
  #output_filename <- "/tmp/output.vcf"
  input_filename <- tempfile("vep-input", fileext = ".vcf")
  #input_filename <- "/tmp/input.vcf"
  output_filename <- tempfile("vep-output", fileext = ".vcf")
  write_tsv(input, file=input_filename)
  print(input)
  vcf_output <- TRUE
  # Note! This will actually produce many rows per variant. These will be joined into one string at the end of this function.
  # The --check-existing will return variants ids such as "rs881144,COSV60975501"
  cmd <- sprintf("~/source/ensembl-vep/vep --vcf --check_existing -i %s -o %s -cache", input_filename, output_filename)
  print(cmd)
  system(cmd)
  result <- read_tsv(output_filename, comment="##")
  #unlink(input_filename)
  #unlink(output_filename)
  #unlink(paste0(output_filename, "_summary.html"))
  
  if (vcf_output) {
    result <- result %>% select("CHR"="#CHROM", POS, INFO)
    names <- rep(NA, 18)
    names[18] <- "rsid"
    result <- result %>% 
      separate(INFO, into=names, remove=FALSE, sep="\\|") %>% 
      mutate(rsid=str_extract(rsid, "rs[0-9]+")) %>%
      rename(Annotation=INFO)
    return(result)
  }
  annotation_variables <- c("Gene", "Feature", "Feature_type", "Consequence", "cDNA_position", "CDS_position", "Protein_position", 
                            "Amino_acids", "Codons", "Extra")
  # #Uploaded_variation     Location        Allele  Gene    Feature Feature_type   >
  # 2_46126027_A/G  2:46126027      G       ENSG00000171132
  result <- result %>%
    separate("#Uploaded_variation", into=c("CHR", "POS", "REF", "ALT"),  sep="[_/]", remove=FALSE) %>%
    #    select(-c(Location, Allele)) %>%
    mutate(across(c(CHR, POS), as.integer)) %>%
    mutate(rsid=str_extract(Existing_variation, "rs[0-9]+")) %>%   # take only the rsid
    unite(col="annotation", all_of(annotation_variables), sep="|") %>% # join the annotation variables together as one string
    group_by(across(-annotation)) %>%
    summarise(annotation=paste0(annotation, collapse="^"))             # one variant can have several annotations, join them all
  print(result)
  return(result)
}

# File: gs://finngen-production-library-green/finngen_R6/finngen_R6_analysis_data/annotations/fin_enriched_genomes_select_columns.txt.gz 
# Column: enrichment_nfsee 
get_enrichment <- function() {
  filename <- "~/data/private/finngen/enrichment/fg_parsed.tsv"
  df <- read_delim(filename, delim=" ", col_names = FALSE)
  colnames(df) <- c("SNPID", "enrichment")
  return(df)
}

get_genotyping_rate <- function() {
  chrs <- c(1:22, "X")
  read_missingness <- function(chr) {
    base <- sprintf("/home/toivoja/FRCBS/blood_health_phewas/results/missingness-chr%s.lmiss", chr)
    filename <- paste0(base, ".lmiss")
    if (!file.exists(filename)) {
      cmd <- sprintf("plink1.9 --bfile /home/toivoja/FRCBS/blood_health_phewas/results/cleaned-deduplicated-chr%s --missing --out %s", chr, base)
      system(cmd)      
    }
    lmiss <- read_delim(filename, delim=" ", trim_ws=TRUE)
    lmiss <- lmiss %>% 
      select(SNPID=SNP, F_MISS) %>%
      mutate(F_MISS=1-F_MISS)
  }
  df <- map_dfr(chrs, read_missingness)
  return(df)
}

get_ld_helper <- function(chr, snips) {
  if (length(snips) <= 1) return(NULL)
  
  if (chr==23) chr <- "X"
  input <- sprintf("~/data/public/1000genomes/supporting/GRCh38_positions/finns.chr%s_GRCh38.genotypes.renamed.reided", chr)
  output <- tempfile("ld")
  snip_file <- tempfile("snips", fileext=".txt")
  write_lines(snips, file=snip_file)
  #cmd <-  sprintf("plink1.9 --vcf %s --vcf-half-call m --r2 --ld-snp-list %s --out %s", input, snip_file, output)
  #  cmd <-  sprintf("plink1.9 --bfile %s --r2 --ld-snp-list %s --out %s", input, snip_file, output)
  # https://groups.google.com/g/plink2-users/c/DO-PEK0lUs8
  cmd <-  sprintf("plink1.9 --bfile %s --r2 inter-chr dprime --ld-window-r2 0.0 --extract %s --out %s", input, snip_file, output)
  
  system(cmd)
  ld_filename <- paste0(output, ".ld")
  if (!file.exists(ld_filename))  {
    cat(sprintf("Error in getting LD for snips in chromosome %s. Snips are:\n", chr))
    print(snips)
    return(NULL)
  }
  # Get rid of the extra spaces
  cleaned_filename <- tempfile("processed", fileext=".txt")
  cmd2 <- sprintf("sed -r 's/^ +//;s/ +$//;s/ +/\\t/g' %s > %s", ld_filename, cleaned_filename)
  system(cmd2)
  
  df <- read_tsv(cleaned_filename)
  df <- df %>% 
    select(SNPID=SNP_A, NEXT_SNPID=SNP_B, R2_TO_NEXT=R2) %>%
    mutate(R2_TO_NEXT=as.character(R2_TO_NEXT))
  return(df)
}

get_ld <- function(df) {
  chrs <- unique(df$CHR)
  lds <- map_dfr(chrs, function(chr) get_ld_helper(chr, df %>% filter(CHR==chr) %>% pull(SNPID)))
  return(lds)
}

add_links <- function(df) {
  format1 <- '=HYPERLINK("https://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=%s", "%s")'
  format2 <- '=HYPERLINK("https://www.snpedia.com/index.php/%s", "%s")'
  format3 <- '=HYPERLINK("https://www.ensembl.org/Homo_sapiens/Variation/Explore?db=core;v=%s", "%s")'
  df <- df %>% mutate(dbSNP   = ifelse(!is.na(rsid), sprintf(format1, rsid, rsid), NA),
                      snipedia= ifelse(!is.na(rsid), sprintf(format2, rsid, rsid), NA),
                      ensembl = ifelse(!is.na(rsid), sprintf(format3, rsid, rsid), NA))
  class(df$dbSNP)    <- c("formula", "character")  # openxlsx::write.xlsx can understand the formula class
  class(df$snipedia) <- c("formula", "character")  # openxlsx::write.xlsx can understand the formula class
  class(df$ensembl)  <- c("formula", "character")  # openxlsx::write.xlsx can understand the formula class
  return(df)
}

# Get the effect sizes of the null model
get_null_model_coefficients <- function(dir, phenotype) {
  filename <- sprintf("%s/results.%s.rda", dir, phenotype)
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
