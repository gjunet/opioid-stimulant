# =============================================================================
# MASTER REPLICATION SCRIPT
#
# Run this file once to execute the complete analysis in the required order.
# Each numbered script is evaluated in a fresh environment so that temporary
# objects from one analysis cannot silently affect a later analysis. Shared
# file locations come from config.R, and shared coding/model functions come
# from R/helpers.R.
# =============================================================================

root <- normalizePath(getwd(), mustWork = TRUE)
if (!file.exists(file.path(root, "config.R"))) stop("Run this file from the publication_code directory.")
source(file.path(root, "config.R"))
# The order matters: scripts 02-07 require the derived file created by script 01,
# and Script 06 verifies the outputs only after all model scripts have run.
scripts <- c("01_build_analysis_data.R", "02_primary_tables.R",
             "03_sensitivity_calendar_period.R", "04_sensitivity_stimulant_class.R",
             "05_sensitivity_code_definition.R",
             "07_sensitivity_alcohol_adjustment.R", "06_verify_outputs.R")
for (script in scripts) {
  message("\n=== Running ", script, " ===")
  # A new child environment retains access to configuration values in the
  # global environment but keeps script-specific objects separate.
  source(file.path(root, "scripts", script), local = new.env(parent = globalenv()))
}
message("\nReplication complete. Tables: ", table_dir)
