load_single <- function(filename) {
  names <- load(filename, verbose=FALSE)
  stopifnot(length(names) == 1)
  return(get(names))
}

ndonor <- function(df) n_distinct(df$donor)

get_finngen <- function() {
  finngen <- load_single("~/proj/interval_prediction/data/preprocessed_progesa.RData")  # This is actually Finngen data combined with progesa
  
  old_count <- nrow(finngen); old_count2 <- ndonor(finngen)
  finngen <- finngen %>% distinct()
  cat(sprintf("FinnGen: Dropped %i / %i donations (%i / %i donors) because of duplicate rows\n",
              old_count - nrow(finngen), old_count, old_count2 - ndonor(finngen), old_count2))
  
  old_count <- nrow(finngen); old_count2 <- ndonor(finngen)
  bad_donors <- finngen %>% count(donor) %>% filter(n>=2) %>% pull(donor)
  finngen <- finngen %>% filter(!(donor %in% bad_donors))
  cat(sprintf("FinnGen: Dropped %i / %i donations (%i / %i donors) that contain contradictory donor information\n",
              old_count - nrow(finngen), old_count, old_count2 - ndonor(finngen), old_count2))
  
  old_count <- nrow(finngen); old_count2 <- ndonor(finngen)
  finngen <- finngen %>%
    select(donor, "height" = "Height..cm.", "weight"="Weight..kg.", "smoking" = "Smoking", gender) %>%
    #rename("height" = "Height..cm.", "weight"="Weight..kg.", "smoking" = "Smoking") %>%
    drop_na() %>%
    #    mutate(smoking = as.character(smoking)) %>%
    mutate(smoking = as.factor(smoking))
  cat(sprintf("FinnGen: Dropped %i / %i donations (%i / %i donors) because of NA values\n",
              old_count - nrow(finngen), old_count, old_count2 - ndonor(finngen), old_count2))
  return(finngen)
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
  "prs_anemia", "PRS Anemia", "numeric", "Anemia PRS",
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
  
