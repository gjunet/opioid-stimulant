# =============================================================================
# SCRIPT 01: BUILD THE ANALYSIS DATASET
#
# Purpose
# Read the ten annual death files, identify deaths involving opioids or
# stimulants, code the demographic adjustment variables and medical-condition
# outcomes, and combine the annual results into one analysis file.
#
# Unit of analysis
# One death certificate. A death is retained when at least one qualifying
# opioid or stimulant code appears as the underlying cause or in any available
# record-axis multiple-cause field.
#
# Privacy
# Source files are read but never modified. The derived RDS remains local and
# contains record-level data; it must not be shared publicly. Logs and tables
# contain aggregate information only.
# =============================================================================

# ---- 1. Load shared functions and confirm that all annual files exist --------
source(file.path(publication_root, "R", "helpers.R"))
log_msg <- new_logger("01_build_analysis_data.log")
required_files <- sprintf(raw_template, years)
if (any(!file.exists(required_files))) stop("Missing input files: ", paste(required_files[!file.exists(required_files)], collapse=", "))

# ---- 2. Define the medical-condition indicators -----------------------------
#
# These regular-expression patterns represent ICD-10 causes associated with
# the prespecified health-condition domains used in the paper. For each retained
# death, the script searches the underlying cause and every available record-axis
# field. A domain is coded 1 when any field begins with one of its listed codes
# and 0 otherwise. The leading ^ means "the code starts here," which prevents a
# match on the same characters appearing later in an unrelated code.
#
# The domains are:
# - cardiac_hypertension: heart disease, selected vascular disease, and
#   hypertension codes listed in the manuscript definition;
# - cancer: malignant neoplasms and the specified HIV-related neoplasm codes;
# - cerebrovascular: cerebrovascular diseases, I60-I69;
# - diabetes_metabolic: diabetes mellitus, E10-E14;
# - hepatobiliary_pancreatic: selected liver, biliary, and pancreatic diseases;
# - renal: selected glomerular, acute/chronic kidney, and other renal diseases;
# - viral_hepatitis: B15-B19;
# - hiv_aids: B20-B24;
# - substance_use_or_mental_health_disorder: F10-F19 diagnoses; and
# - accidental_injury_other_than_drug_overdose: external injury codes after
#   excluding drug-overdose intent codes and later excluding suicide/homicide.
#
# These flags describe codes recorded on the death certificate. They should not
# be interpreted as a complete clinical history or as conditions caused by the
# opioid/stimulant exposure.
condition_patterns <- c(
  cardiac_hypertension = "^(I0[0-9]|I1[0-3]|I20|I21|I22|I24|I25|I30|I31|I33|I40|I50|I26|I27|I28|I34|I35|I36|I37|I38|I42|I43|I44|I45|I46|I47|I48|I49|I51|I70|I71|I10|I11|I12|I15)",
  cancer = "^(B212|B218|C0[0-9]|C1[0-9]|C2[0-9]|C3[0-9]|C4[0-9]|C5[0-9]|C6[0-9]|C7[0-9]|C8[0-9]|C9[0-7])",
  cerebrovascular = "^I6[0-9]", diabetes_metabolic = "^(E10|E11|E12|E13|E14)",
  hepatobiliary_pancreatic = "^(K70|K73|K74|K76|K72|K80|K81|K82|K85|K86)",
  renal = "^(N00|N01|N02|N03|N04|N05|N06|N07|N17|N18|N19|N25|N26|N27|N28)",
  viral_hepatitis = "^(B15|B16|B17|B18|B19)", hiv_aids = "^(B20|B21|B22|B23|B24)",
  substance_use_or_mental_health_disorder = "^F1",
  accidental_injury_other_than_drug_overdose = "^(W|V|Y0[0-9]|Y1[5-9]|Y[2-9]|X1|X2|X3|X6|X7|X8|X9|X0)")
# Search both the underlying cause and up to 20 multiple-cause record-axis fields.
icd_fields <- c("ucod_icd10", paste0("recaxis", 1:20))

# ---- 3. Process one annual mortality file at a time --------------------------
# Processing by year limits memory use. Each annual result is saved temporarily
# in the parts list and combined only after every year has been coded.
parts <- list()
for (year in years) {
  # Read one authorized annual file and identify the cause fields it contains.
  raw <- readRDS(sprintf(raw_template, year))
  fields <- intersect(icd_fields, names(raw))
  if (!"ucod_icd10" %in% fields) stop(year, ": ucod_icd10 is required")
  # Remove punctuation and capitalization differences from all cause codes,
  # then concatenate the fields for the drug-involvement search. Spaces remain
  # between fields so no synthetic code can be created at a field boundary.
  normalized <- lapply(raw[fields], normalize_icd)
  combined <- do.call(paste, c(normalized, sep=" "))

  # Inclusive drug definition used in the primary paper:
  # opioid = F11, T40.0-T40.4, T40.6, or R78.1;
  # cocaine = F14, T40.5, or R78.2; and
  # other psychostimulant = F15 or T43.6.
  # T40.5 is cocaine and is therefore stimulant, not opioid, involvement.
  # Maternal-exposure P04 codes are intentionally absent.
  opioid <- grepl("F11|T40[0-46]|R781", combined)
  cocaine <- grepl("F14|T405|R782", combined)
  psychostimulant <- grepl("F15|T436", combined)
  stimulant <- cocaine | psychostimulant
  # Keep deaths involving at least one of the two broad drug classes.
  keep <- opioid | stimulant

  # Restricted sensitivity definition. It removes all F11/F14/F15 mental and
  # behavioral disorder codes and retains poisoning or blood-finding codes only.
  # Some deaths retained by the inclusive definition have no restricted code;
  # their restricted_group is set to missing and Sensitivity 3 excludes them.
  restricted_opioid <- grepl("T40[0-46]|R781", combined[keep])
  restricted_stimulant <- grepl("T405|T436|R782", combined[keep])
  # Search each cause field separately for every medical-condition domain. The
  # logical OR means that one matching field is sufficient to code the domain 1.
  norm_kept <- lapply(normalized, function(x) x[keep])
  flags <- sapply(names(condition_patterns), function(name) {
    out <- rep(FALSE, sum(keep))
    for (field in norm_kept) out <- out | grepl(condition_patterns[name], field)
    as.integer(out)
  })
  flags <- as.data.frame(flags); names(flags) <- names(condition_patterns)
  # The broad external-injury pattern can capture intentional injuries. Set the
  # accidental-injury indicator to 0 when manner of death is suicide or homicide.
  manner <- as.character(get_col(raw,"mannerdeath"))[keep]
  flags$accidental_injury_other_than_drug_overdose[manner %in% c("2","3")] <- 0L
  # Record whether any cause field contains a prespecified alcohol-attributable
  # code. This variable is descriptive and is not an adjustment covariate.
  alcohol <- rep(FALSE,sum(keep)); for(field in norm_kept) alcohol <- alcohol | grepl("^(F10|K70|X45|T51)",field)

  # Assemble the annual analytic extract. comparison_group uses the inclusive
  # drug definition; restricted_group uses the poisoning/blood-only definition.
  # Separate cocaine and psychostimulant flags support Sensitivity 2.
  y <- rep(year,sum(keep))
  part <- data.frame(year=y, comparison_group=drug_group(opioid[keep],stimulant[keep]),
    restricted_group=drug_group(restricted_opioid,restricted_stimulant),
    opioid=opioid[keep],cocaine=cocaine[keep],psychostimulant=psychostimulant[keep],
    age_years=decode_age(get_col(raw,"age")[keep]), sex=decode_sex(get_col(raw,"sex")[keep]),
    race4=decode_race(y,get_col(raw,"race")[keep],get_col(raw,"racerec5")[keep],get_col(raw,"racerec40")[keep]),
    hisp3=decode_hispanic(y,get_col(raw,"hispanic")[keep],get_col(raw,"hisp_origin")[keep]),
    alcohol_any=as.integer(alcohol),stringsAsFactors=FALSE)
  # Add the medical-condition indicators and retain the annual extract.
  part <- cbind(part,flags); parts[[as.character(year)]] <- part
  log_msg(year, ": source N=",format(nrow(raw),big.mark=","),"; inclusive N=",format(sum(keep),big.mark=","),"; restricted N=",format(sum(restricted_opioid|restricted_stimulant),big.mark=","))
  rm(raw,normalized,combined,norm_kept,flags,part); gc(verbose=FALSE)
}
# ---- 4. Combine years and construct the four manuscript outcomes -------------
d <- do.call(rbind,parts); rownames(d) <- NULL

# Cardiovascular reproduces the cardiac/hypertension domain. "Other medical"
# is 1 when any cancer, diabetes/metabolic, liver/biliary/pancreatic, renal,
# viral-hepatitis, or HIV/AIDS code is present. "No additional" is 1 only when
# cardiovascular, cerebrovascular, and other medical are all absent.
d$cardiovascular <- d$cardiac_hypertension
d$other_medical <- as.integer(d$cancer | d$diabetes_metabolic | d$hepatobiliary_pancreatic | d$renal | d$viral_hepatitis | d$hiv_aids)
d$no_additional <- as.integer(!(d$cardiovascular | d$cerebrovascular | d$other_medical))
# ---- 5. Save the local derived file and a non-record-level variable inventory -
# The RDS is needed by scripts 02-06 but is excluded by .gitignore. The codebook
# lists variable names and storage classes without exposing death-level values.
saveRDS(d,analytic_path)
codebook <- data.frame(variable=names(d),class=vapply(d,function(x)class(x)[1],character(1)),stringsAsFactors=FALSE)
write.csv(codebook,file.path(derived_dir,"analytic_variable_codebook.csv"),row.names=FALSE)
log_msg("Saved derived analytic data; N=",format(nrow(d),big.mark=","))
