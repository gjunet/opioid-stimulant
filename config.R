# =============================================================================
# CONFIGURATION
#
# This file tells every analysis script where to find the source mortality
# files and where to write derived data, tables, diagnostics, and logs. Most
# users only need to edit this file if their data are stored outside the
# publication_code folder. Environment variables take precedence over the
# default subfolders, which is useful on secure servers where record-level
# mortality data must remain in a protected location.
# =============================================================================

# Treat the directory from which run_all.R is launched as the replication root.
publication_root <- normalizePath(getwd(), mustWork = TRUE)

# Input, working-data, and output locations. Script 01 reads raw_dir and writes
# a local analysis file to derived_dir. The remaining scripts read that local
# analysis file and write aggregate results only.
raw_dir <- Sys.getenv("STIM_OPIOID_RAW_DIR", file.path(publication_root, "data", "raw"))
derived_dir <- Sys.getenv("STIM_OPIOID_DERIVED_DIR", file.path(publication_root, "data", "derived"))
output_dir <- Sys.getenv("STIM_OPIOID_OUTPUT_DIR", file.path(publication_root, "output"))
table_dir <- file.path(output_dir, "tables")
diagnostic_dir <- file.path(output_dir, "diagnostics")
log_dir <- file.path(publication_root, "logs")
# Fixed manuscript study period and required annual filename pattern.
years <- 2015:2024
raw_template <- file.path(raw_dir, "processed_mort%d.rds")
analytic_path <- file.path(derived_dir, "stim_opioid_analysis_2015_2024.rds")
# Create output folders if they do not already exist. This does not change or
# write to the source mortality files.
dir.create(derived_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(diagnostic_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
