# =============================================================================
# SCRIPT 06: VERIFY THE REPLICATION
#
# This final script does not fit substantive models. It confirms that every
# promised table exists and every fitted model reported convergence. It also
# records the R session and creates a local file manifest. The public repository
# contains no expected-results file or generated numerical output.
# =============================================================================

# ---- 1. Confirm that every primary and sensitivity table was created ---------
source(file.path(publication_root,"R","helpers.R")); log_msg <- new_logger("06_verify_outputs.log")
expected_files <- c("table1_chang_2015_2024.csv","S1_descriptive_2015_2024.csv","S2_unadjusted_2015_2024.csv","S3_adjusted_2015_2024.csv","S4_pairwise_cardiovascular.csv","S5_pairwise_cerebrovascular.csv","S6_pairwise_other_medical.csv","S7_pairwise_no_additional.csv","eTable_icd10_drug_codes.csv","sensitivity1_period_coefficients.csv","sensitivity1_period_change_tests.csv","sensitivity2_stimulant_class_samples.csv","sensitivity2_stimulant_class_models.csv","sensitivity3_inclusive_restricted_models.csv","sensitivity3_sample_sizes.csv","sensitivity3_outcome_prevalence.csv")
missing <- expected_files[!file.exists(file.path(table_dir,expected_files))]
if(length(missing)) stop("Missing expected outputs: ",paste(missing,collapse=", "))
# ---- 2. Read all diagnostic files and stop after a convergence failure -------
# Stopping here prevents a nominally complete replication from being reported
# when one or more logistic models failed to reach a solution.
diag_files <- list.files(diagnostic_dir,pattern="diagnostics\\.csv$",full.names=TRUE)
diag_data <- do.call(rbind,lapply(diag_files,read.csv,stringsAsFactors=FALSE))
if(any(!diag_data$converged)) stop("At least one model did not converge. See output/diagnostics.")
# ---- 3. Write an output manifest and software-version record -----------------
# These files are generated locally and ignored by Git.
manifest <- data.frame(file=expected_files,path=file.path("output","tables",expected_files),bytes=file.info(file.path(table_dir,expected_files))$size,stringsAsFactors=FALSE)
write.csv(manifest,file.path(output_dir,"output_manifest.csv"),row.names=FALSE)
session <- capture.output(sessionInfo()); writeLines(session,file.path(output_dir,"sessionInfo.txt"))
log_msg("All expected files present and all models converged.")
