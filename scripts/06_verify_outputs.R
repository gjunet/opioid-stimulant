# =============================================================================
# SCRIPT 06: VERIFY THE REPLICATION
#
# This final script does not fit substantive models. It confirms that every
# promised table exists and every fitted model reported convergence. It also
# records the R session and creates a local file manifest. The public repository
# contains no expected-results file or generated numerical output.
# =============================================================================

# ---- 1. Confirm that every primary, sensitivity, and figure file was created -
source(file.path(publication_root,"R","helpers.R")); log_msg <- new_logger("06_verify_outputs.log")
expected_tables <- c("table1_chang_2015_2024.csv","S1_descriptive_2015_2024.csv","S2_unadjusted_2015_2024.csv","S3_adjusted_2015_2024.csv","S4_pairwise_cardiovascular.csv","S5_pairwise_cerebrovascular.csv","S6_pairwise_other_medical.csv","S7_pairwise_no_additional.csv","eTable_icd10_drug_codes.csv","sensitivity1_period_coefficients.csv","sensitivity1_period_change_tests.csv","sensitivity2_stimulant_class_samples.csv","sensitivity2_stimulant_class_models.csv","sensitivity3_inclusive_restricted_models.csv","sensitivity3_sample_sizes.csv","sensitivity3_outcome_prevalence.csv","sensitivity4_alcohol_adjustment_models.csv","Figure1_plotted_estimates.csv","Figure1_plotted_estimates_opioid_reference.csv")
expected_figures <- c("Figure1_adjusted_odds_ratios.png","Figure1_adjusted_odds_ratios.svg","Figure1_adjusted_odds_ratios_linear_scale.png","Figure1_adjusted_odds_ratios_linear_scale.svg","Figure1_adjusted_logit_coefficients.png","Figure1_adjusted_logit_coefficients.svg","Figure1_adjusted_odds_ratios_opioid_reference.png","Figure1_adjusted_odds_ratios_opioid_reference.svg","Figure1_adjusted_odds_ratios_linear_scale_opioid_reference.png","Figure1_adjusted_odds_ratios_linear_scale_opioid_reference.svg","Figure1_adjusted_logit_coefficients_opioid_reference.png","Figure1_adjusted_logit_coefficients_opioid_reference.svg")
missing_tables <- expected_tables[!file.exists(file.path(table_dir,expected_tables))]
missing_figures <- expected_figures[!file.exists(file.path(figure_dir,expected_figures))]
if(length(missing_tables) || length(missing_figures)) stop(
  "Missing expected outputs: ",paste(c(missing_tables,missing_figures),collapse=", "))
# ---- 2. Read all diagnostic files and stop after a convergence failure -------
# Stopping here prevents a nominally complete replication from being reported
# when one or more logistic models failed to reach a solution.
diag_files <- list.files(diagnostic_dir,pattern="diagnostics\\.csv$",full.names=TRUE)
diagnostics_converged <- vapply(diag_files,function(path) {
  diagnostic <- read.csv(path,stringsAsFactors=FALSE)
  if(!"converged" %in% names(diagnostic)) stop("Missing convergence field in ",basename(path))
  all(diagnostic$converged)
},logical(1))
if(any(!diagnostics_converged)) stop("At least one model did not converge. See output/diagnostics.")
# ---- 3. Write an output manifest and software-version record -----------------
# These files are generated locally and ignored by Git.
manifest <- rbind(
  data.frame(file=expected_tables,path=file.path("output","tables",expected_tables),bytes=file.info(file.path(table_dir,expected_tables))$size,stringsAsFactors=FALSE),
  data.frame(file=expected_figures,path=file.path("output","figures",expected_figures),bytes=file.info(file.path(figure_dir,expected_figures))$size,stringsAsFactors=FALSE)
)
write.csv(manifest,file.path(output_dir,"output_manifest.csv"),row.names=FALSE)
session <- capture.output(sessionInfo()); writeLines(session,file.path(output_dir,"sessionInfo.txt"))
log_msg("All expected files present and all models converged.")
