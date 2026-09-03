# Stimulant–Opioid Comorbidity Study: Replication Code

This folder contains the R code needed to recreate the manuscript tables and four post hoc sensitivity analyses for deaths occurring during 2015–2024.

The mortality microdata are not distributed with this repository. Researchers must obtain authorized access to the National Center for Health Statistics Multiple Cause of Death data and prepare the annual input files described below.

## Software

- R 4.2 or later
- `ggplot2` is required to reproduce Figure 1; the remaining analyses use base R.

## Input files

Place one RDS file per year in `data/raw/`:

```text
data/raw/processed_mort2015.rds
...
data/raw/processed_mort2024.rds
```

Each file must contain one death per row and these variables:

- `ucod_icd10`: underlying cause of death;
- `recaxis1` through `recaxis20`: record-axis multiple-cause fields (fewer fields are accepted, but `ucod_icd10` is required);
- `age`, `sex`, `race`, `racerec5`, `racerec40`, `hispanic`, `hisp_origin`, and `mannerdeath`.

Variable names must be lowercase. ICD-10 values may contain decimal points or spaces; the code normalizes them before matching. Missing demographic variables may be present as all-missing columns when a coding-era field does not apply.

The annual RDS files are treated as read-only. No record-level data are written by the table scripts. Script 01 creates a derived analytic RDS locally so later scripts do not repeatedly scan the source files; this derived file must not be uploaded to a public repository.

## Running the replication

From this folder, run:

```r
source("run_all.R")
```

To use a different input or output location, edit `config.R` or set these environment variables before running R:

```text
STIM_OPIOID_RAW_DIR
STIM_OPIOID_DERIVED_DIR
STIM_OPIOID_OUTPUT_DIR
```

Tables are written to `output/tables/`, figures to `output/figures/`, diagnostics to `output/diagnostics/`, and timestamped logs to `logs/`.

## Rendering the data dictionary as PDF

Pandoc and a LaTeX installation with XeLaTeX are required. From the repository root, run:

```sh
pandoc DATA_DICTIONARY.md \
  --include-in-header=pandoc-header.tex \
  --pdf-engine=xelatex \
  -o DATA_DICTIONARY.pdf
```

The condition table uses short display labels so its first column wraps correctly in LaTeX. Exact R variable names appear in a full-width list immediately after the table.

## Analysis sequence

1. `scripts/01_build_analysis_data.R` codes the inclusive and restricted drug definitions, the separate cocaine and other-psychostimulant indicators, demographics, and outcomes.
2. `scripts/02_primary_tables.R` creates manuscript Table 1, Appendix Tables S1–S7, and the ICD-10 code table.
3. `scripts/03_sensitivity_calendar_period.R` compares 2015–2019 with 2020–2024 using period-specific models and formal interaction tests.
4. `scripts/04_sensitivity_stimulant_class.R` analyzes cocaine and other psychostimulants separately.
5. `scripts/05_sensitivity_code_definition.R` compares the inclusive and restricted drug-involvement definitions.
6. `scripts/06_verify_outputs.R` checks that expected files were generated and that every model converged, then writes a local replication manifest.
7. `scripts/07_sensitivity_alcohol_adjustment.R` compares the primary models with otherwise identical models that add death-certificate alcohol involvement as a covariate.
8. `scripts/08_figure1_forest_plots.R` creates Figure 1 on both the adjusted odds-ratio and logit-coefficient scales. It also creates companion versions with opioid-only deaths as the reference group.

## Figure 1 outputs

The primary Figure 1 files use stimulant-only deaths as the reference group and are saved as `Figure1_adjusted_odds_ratios.*` and `Figure1_adjusted_logit_coefficients.*`. Companion files ending in `_opioid_reference` reproduce the alternate opioid-only reference specification. Each point is labeled with its estimate and 95% confidence interval. Reference-group and adjustment information appears beneath the plotting area rather than as a title or subtitle.

## How to read the scripts

Each numbered script begins with a plain-language statement of its purpose and then uses numbered comment headings. Read the paragraph under a heading before reading the R statements that follow it. The comments explain the analytic reason for the step, the records affected, the reference category, and how the resulting quantity should be interpreted.

`R/helpers.R` contains operations shared across analyses. For example, it standardizes ICD-10 text, translates NCHS demographic variables that changed across years, assigns the three mutually exclusive drug groups, and fits the common adjusted logistic model. Centralizing these functions prevents the sensitivity analyses from silently using different demographic or modeling rules.

`DATA_DICTIONARY.md` defines every derived exposure, demographic variable, condition indicator, and modeled outcome in plain language.

## Supplementary appendix

The complete publication supplement is available in two formats:

- [`docs/Supplementary_Appendix.md`](docs/Supplementary_Appendix.md), an accessible and version-controlled Markdown source;
- [`docs/Supplementary_Appendix.pdf`](docs/Supplementary_Appendix.pdf), the formatted publication version.

The supplement contains the data dictionary and all four post hoc sensitivity analyses: calendar period, separate cocaine and other-psychostimulant groups, inclusive versus restricted drug-code definitions, and adjustment for alcohol involvement.

The cause-code searches use regular expressions. After punctuation is removed, `^I6[0-9]`, for example, means that a code starts with I6 and ends that three-character stem with any digit from 0 through 9; this identifies I60-I69 cerebrovascular codes. The comments above `condition_patterns` in Script 01 explain every health-condition domain. A condition flag means that at least one qualifying code was recorded as the underlying or a multiple cause of death. It does not establish onset, severity, or causation by the drug exposure.

The code uses short internal names such as `x`, `dd`, and `fit` only within a small analysis block. In those blocks, `x` or `dd` is the current analysis dataset and `fit` is the fitted logistic regression. Longer-lived objects use descriptive names such as `restricted_group`, `psychostimulant`, and `period`.

## Drug definitions

The inclusive opioid definition is F11, T40.0–T40.4, T40.6, or R78.1. The inclusive stimulant definition is F14, F15, T40.5, T43.6, or R78.2. T40.5 is classified as stimulant involvement and is not counted as opioid involvement.

The restricted definition retains only T40.0–T40.6, T43.6, R78.1, and R78.2, with T40.5 assigned to stimulant involvement. It excludes all F11, F14, and F15 mental and behavioral disorder codes.

Maternal-exposure codes P04.14, P04.16, and P04.41 are excluded from both definitions.

## Statistical models

The four binary outcomes are cardiovascular disease, cerebrovascular disease, other medical conditions, and no additional condition. Logistic regression models adjust for age in 10-year units, sex, race, and Hispanic origin. The primary reference group is stimulant-only deaths. Missing model variables are handled by complete-case analysis through R's `glm` default.

All sensitivity analyses were specified after the primary analysis and should be described as post hoc.

## Reproducibility notes

Results can differ if the source files use different NCHS revisions, race/Hispanic-origin recodes, cause-field representations, or record inclusion rules. The verification script stops when required outputs are missing or a model does not converge. Generated tables, diagnostics, logs, and derived files are excluded from the public repository; they are created locally when an authorized user runs the code.

## Suggested citation

When publishing the repository, replace the placeholder citation information in `CITATION.cff` with the final article citation and repository DOI.
