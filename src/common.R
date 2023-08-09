
library(tidyverse)

load_single <- function(filename) {
  names <- load(filename, verbose=FALSE)
  stopifnot(length(names) == 1)
  return(get(names))
}

ndonor <- function(df) n_distinct(df$donor)

get_finngen <- function(na.rm = TRUE) {
  finngen <- load_single("~/proj/interval_prediction/data/preprocessed_progesa.RData")  # This is actually Finngen data combined with progesa

  # Concentrate only on relevant variables  
  finngen <- finngen %>%
    select(donor, "height" = "Height..cm.", "weight"="Weight..kg.", "smoking" = "Smoking", gender) %>%
    #rename("height" = "Height..cm.", "weight"="Weight..kg.", "smoking" = "Smoking") %>%
    #    mutate(smoking = as.character(smoking)) %>%
    mutate(smoking = as.factor(smoking))

  # Drop duplicate rows
  old_count <- nrow(finngen); old_count2 <- ndonor(finngen)
  finngen <- finngen %>% distinct()
  cat(sprintf("FinnGen: Dropped %i / %i donations (%i / %i donors) because of duplicate rows\n",
              old_count - nrow(finngen), old_count, old_count2 - ndonor(finngen), old_count2))
  
  # Drop donors that have conflicting donor information
  old_count <- nrow(finngen); old_count2 <- ndonor(finngen)
  bad_donors <- finngen %>% count(donor) %>% filter(n>=2) %>% pull(donor)
  finngen <- finngen %>% filter(!(donor %in% bad_donors))
  cat(sprintf("FinnGen: Dropped %i / %i donations (%i / %i donors) that contain contradictory donor information\n",
              old_count - nrow(finngen), old_count, old_count2 - ndonor(finngen), old_count2))
  
  
  if (na.rm) {
    old_count <- nrow(finngen); old_count2 <- ndonor(finngen)
    finngen <- finngen %>%
      drop_na() 
    cat(sprintf("FinnGen: Dropped %i / %i donations (%i / %i donors) because of NA values\n",
                old_count - nrow(finngen), old_count, old_count2 - ndonor(finngen), old_count2))
  }
  
  return(finngen)
}

# Full finngen includes also FinDonor
get_full_finngen <- function() {
  finngen2 <- get_finngen(na.rm=FALSE)
  finngen2 <- finngen %>% mutate(smoking = fct_explicit_na(smoking, na_level = NA))
  fd_filename <- "~/proj/interval_prediction/data/DataFromMuriel/20181129/r02ds.donorData.rdata"
  findonor <- load_single(fd_filename)
  findonor <- findonor %>% anti_join(finngen2, by="donor") # There are four donors that are in both datasets.
  # In that case, use the information in FinnGen
  findonor2 <- findonor %>% 
    select(donor, weight, height, smoking=QR54, gender=Gender) %>%
    mutate(smoking = fct_recode(smoking, 
                                `yes, occassionally` = "sometimes",
                                `yes, regularly` = "daily"),
           smoking = factor(smoking, ordered = FALSE),
           smoking = fct_explicit_na(smoking, na_level = NA))
  finngen <- bind_rows(finngen2, findonor2)  
}

donation_descript <- tibble(Variable = c("donor", "Hb", "days_to_previous_fb", "age", "previous_Hb_def", 
                                "year", "warm_season", "consecutive_deferrals", "recent_donations",
                                "recent_deferrals", "hour", 
                                "previous_Hb", "Hb_first", "Hb_deferral", "sex"), 
                   Pretty = c("Donor ID", "Hemoglobin", "Days to previous full blood donation", "Age", "Previous event was deferral", 
                              "Year", "Warm season", "Consecutive deferrals", "Donations in last two years", 
                              "Deferrals in last two years", "Hour", 
                              "Previous Hb", "First Hb", "Deferral", "Sex"),
                   Type = c("factor", "numeric", "numeric (int)", "numeric", "boolean",
                            "numeric (int)", "boolean", "numeric (int)", "numeric (int)", "numeric (int)", "numeric",
                            "numeric", "numeric", "boolean", "factor"),
                   Explanation = c("Donor identifier",
                                   "Amount of Hemoglobin",
                                   "Time (in days) between Hb measurement and previous full blood donation event",
                                   "Age of donor",
                                   "Indicates whether the donor had low hemoglobin at previous donation event",
                                   "Year of donation",
                                   "True if donation was given in April-September",
                                   "Number of times the donor has been deferred due to low hemoglobin since last succesful whole blood donation",
                                   "Number of donations in the last two years",
                                   "Number of deferrals in the last two years",
                                   "Time of day when donation was given as hours (e.g. 13:45 = 13.75)",
                                   "Hb value at previous measurement (dynamic linear mixed model)",
                                   "Hb value at first donation of this donor (linear mixed model)",
                                   "Hemoglobin below deferral threshold",
                                   "Sex of the donor")
)


donor_descript <- tibble(
  Variable    = c("smoking", "height", "weight", "RNF43_mutant", "prs", "FERRITIN_FIRST", "FERRITIN_LAST", "one_deferral", "label"),
  Pretty      = c("Smoking", "Height", "Weight", "RNF43", "Polygenic score", "First ferritin", "Last ferritin", "At least one deferral", "Partition label"),
  Type        = c("boolean", "numeric", "numeric", "boolean", "numeric", "numeric", "numeric", "numeric (int)", "factor"),
  Explanation = c("Does the person smoke", "Height of the donor", "Weight of the donor", 
                  "Mutation at RNF43 gene in chromosome 17 position 58358769", "Polygenic risk score for hemoglobin", 
                  "First measured ferritin value", "Last measured ferritin value", "At least one deferral",
                  "The donors are partitioned into train, validate, and test sets")
)

extra_descript <- tribble(
  ~Variable, ~Pretty, ~Type, ~Explanation,
  "prs_anemia", "PRS IDA", "numeric", "Iron deficiency anemia PRS",
  "prs_ferritin", "PRS Ferritin", "numeric", "Ferritin PRS",
  "prs_hemoglobin", "PRS Hemoglobin", "numeric", "Hemoglobin PRS",
  "snp_17_58358769", "SNP 17:58358769", "numeric", "SNP 17:58358769",  
  "snp_6_32617727", "SNP 6:32617727", "numeric", "SNP 6:32617727",
  "snp_15_45095352", "SNP 15:45095352", "numeric", "SNP 15:45095352",
  "snp_1_169549811", "SNP 1:169549811", "numeric", "SNP 1:169549811",
  "female", "Female", "boolean", "Is the donor female",
  "nb_donat", "Number of lifetime donations", "numeric (int)", "Number of lifetime donations",
  "bmi", "BMI", "numeric", "Body Mass Index",
  "blood_donor", "Is blood donor", "boolean", "Is the individual a blood donor"
)

descript <- bind_rows(donation_descript, donor_descript, extra_descript)

# Unfinished
to_pretty <- function(df, descript) {  
  df %>%
    left_join(descript, by=c("Variable")) %>%
    mutate(Pretty = ifelse(is.na(Pretty), Variable, Pretty)) %>%
    mutate(Pretty = as.factor(Pretty))
}

# Converts column names of a dataframe 'df' to pretty format using dataframe 'description'.
pretty_col_names <- function(df, description) {
  conversion <- description %>% #filter(Variable %in% all_of(old_names)) %>% 
    select(Pretty, Variable) %>% deframe()
  df %>% rename(any_of(conversion))
}

#to_pretty_vector <- function(v, descript) {
#  df <- tibble(Variable = v)
#  to_pretty(df, descript)$Pretty
#}

# Make elements of vector pretty
to_pretty_vector <- function(v, description) {
  conversion <- description %>%
    select(Variable, Pretty) %>% deframe()
  recode(v, !!!conversion)
}

# Constructor of the MyLogger class
# S3 class
new_logger <- function(prefix="", file="", silent=FALSE) {
  stopifnot(is.character(prefix))
  structure(list(prefix=prefix, file=file, silent=silent), class = "MyLogger")
}

print.MyLogger <- function(object, msg) {
  if (! object$silent) {
    cat(object$prefix, msg, file=object$file, append=TRUE)
  }
}


closest_snip <- function(df, pattern, radius=NULL) {
  stopifnot(nrow(pattern) == 1)
  result <- df %>% filter(chr==pattern$chr) %>%
    mutate(dist = abs(bp - pattern$bp))
  if (is.null(radius)) {
    result <- result %>% slice_min(order_by = dist, n=1, with_ties = FALSE)   # return closest snip
  } else {
    result <- result %>% filter(dist <= radius, p < 5e-8)   # return all significant snips within radius
  }
  result %>% relocate(chr)
}
get_peaks <- function(df, centers) {
  helper <-  function(...) {
    center <- tibble(...)
    result <- closest_snip(df, center, radius=radius) %>% mutate(center_chr=center$chr, center_bp=center$bp)
    return(result)
  }
  pmap_dfr(centers, helper)
}
get_closest_snips <- function(df, centers) {
  helper <-  function(...) {
    center <- tibble(...)
    result <- closest_snip(df, center) %>% mutate(orig_chr=center$chr, orig_bp=center$bp)
    return(result)
  }
  pmap_dfr(centers, helper)
}
jaccard_index <- function(set1, set2) {
  i <- length(intersect(set1, set2))
  u <- length(union(set1, set2))
  return(list(numerator=i, denominator=u, fraction=i/u))
}
inclusion <- function(set1, set2) {
  i <- length(intersect(set1, set2))
  u <- length(set1)
  return(list(numerator=i, denominator=u, fraction=i/u))
}

find_lead_snips <- function(tmp, radius=1e6) {
  #tmp <- L2$deferral
  # Slide a 2*size bp window over the dataframe rows. Returns a dataframe that
  # contains for each window a row with minimum p value
  find_lead_snips_helper <- function(df, key) {
    #size <- 1e6
    lead_snips <- slide_index_dfr(df, df$bp, ~slice_min(.x, order_by = p, n=1, with_ties = FALSE), 
                                  .before=radius, .after=radius)
    lead_snips <- lead_snips %>% distinct()
    lead_snips$chr <- key$chr
    return(lead_snips)
  }
  find_lead_snips_better <- function(df, key) {
    chr <- key$chr
    #print(chr)
    min_p <- df[1, "p"]
    min_pos <- df[1, "bp"]
    reported <- FALSE
    result <- c()
    for (i in 1:nrow(df)) {
      current_pos <- df[i, "bp"]
      current_p <- df[i, "p"]
      if (current_pos - min_pos > radius) {
        result <- c(result, min_pos)
        #reported <- TRUE
        min_p <- current_p
        min_pos <- current_pos
      } else if (current_p < min_p) {
        min_p <- current_p
        min_pos <- current_pos
        #reported <- FALSE
      }
    }
    result <- c(result, min_pos)
    #print(result)
    lead_snips <- df %>% filter(bp %in% result) %>% mutate(chr=chr)
    return(lead_snips)    
  }
  find_lead_snips_best <- function(df, key) {
    df$minp <- slide_index_dbl(df, df$bp, ~min(.x$p), .before=radius, .after=radius)
    df$chr <- key$chr
    return(df %>% filter(p == minp))
  }
  #tmp2 <- tmp %>% filter(chr==1)
  #df <- find_lead_snips(tmp2)
  res <- tmp %>% filter(p < 5e-8) %>% group_by(chr) %>% group_map(find_lead_snips_best) %>% bind_rows() %>%
    relocate(chr)
  res
}

get_gene_names <- function(all_gwas) {
  use_vep <- TRUE
  significant <- all_gwas %>% 
    filter(p < 5e-8)
  if (use_vep) {
    if (recompute_annotation) {
      header <- str_split("Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|SYMBOL_SOURCE|HGNC_ID|CLIN_SIG|SOMATIC|PHENO", "\\|")[[1]]
      snips_to_be_identified <- significant %>%
        #head(100) %>%   # for testing purposes
        separate(snp, c("dummy1","dummy2","Allele1","Allele2"), sep="[_:]", remove=FALSE) %>% 
        select(CHR=chr, POS=bp, ID=snp, Allele1, Allele2)
      genes2 <- get_annotations(snips_to_be_identified)
      genes3 <- genes2 %>% 
        mutate(Annotation = str_remove(Annotation, "^CSQ=")) %>% 
        separate_rows(Annotation, sep=",") %>% 
        separate(Annotation, 
                 #sprintf("X%i", 1:26),
                 header,
                 sep="\\|")
      saveRDS(genes3, printf("%s/parsed_annotation.rds", result_base))
    } else {
      genes3 <- readRDS(sprintf("%s/parsed_annotation.rds", result_base))
    }
    genes <- genes3 %>% 
      select(chr=CHR, bp=POS, gene=SYMBOL) %>% 
      distinct() %>%
      filter(gene != "")
    
    ################################
    # Add some missing names by hand
    ################################
    
    # These are without names:
    x1 <- significant %>% filter(phenotype == "hemoglobin") %>% filter(chr==4) %>% slice_min(order_by=p, n=1)
    x2 <- significant %>% filter(phenotype == "ferritin")   %>% filter(chr==6) %>% slice_min(order_by=p, n=1)
    x3 <- significant %>% filter(phenotype == "ferritin")   %>% filter(chr==1) %>% slice_min(order_by=p, n=1)
    # Use the names of the closest snips for them
    y1 <- closest_snip(genes, x1)
    y2 <- closest_snip(genes, x2)
    y3 <- closest_snip(genes, x3)
    genes <- genes %>% 
      add_row(chr=x1$chr, bp=x1$bp, gene=y1$gene) %>%
      add_row(chr=x2$chr, bp=x2$bp, gene=y2$gene) %>%
      add_row(chr=x3$chr, bp=x3$bp, gene=y3$gene)
    
    genes_unique <- genes %>%
      group_by(chr, bp) %>% 
      slice_head(n=1) %>%  # Just select the first gene name for the position. This is quite random.
      ungroup()
  } else {  # use plink
    snips_to_be_identified <- significant %>% select(chr, bp) %>% distinct()
    genes <- get_genes_plink(snips_to_be_identified)
    genes <- genes %>% select(-size)
  }
  genes
}

add_gene_names <- function(df, genes) {
  df %>% left_join(genes, by=c("chr", "bp")) %>%
    replace_na(list(gene=""))
  #    replace_na(list(gene="<noname>"))
}