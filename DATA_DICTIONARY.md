# Derived Analysis Data Dictionary

Script 01 creates a local death record-level RDS datafile with one row per death for all cases meeting the inclusive opioid/stimulant definition. Note that the data is not included with these code but can be obtained from [CDC] (https://www.cdc.gov/comec/data-systems/nvss-mortality-data.html).  
Coding of the raw death data can be found at [NEBR]. ([https://www.nber.org/research/data/mortality-data-vital-statistics-nchs-multiple-cause-death-data](https://www.nber.org/research/data/mortality-data-vital-statistics-nchs-multiple-cause-death-data))

## Identification and exposure variables

| Variable | Meaning |
|---|---|
| `year` | Calendar year of death, 2015–2024. |
| `comparison_group` | Three mutually exclusive groups under the inclusive definition: `stimulant_no_opiate`, `opiate_stimulant`, or `opiate_no_stimulant`. |
| `restricted_group` | The same three groups reconstructed from poisoning and toxicologic blood-finding codes only. Missing means the death had an inclusive F11/F14/F15 code but no qualifying restricted code. |
| `opioid` | At least one F11, T40.0–T40.4, T40.6, or R78.1 code was recorded. |
| `cocaine` | At least one F14, T40.5, or R78.2 code was recorded. |
| `psychostimulant` | At least one F15 or T43.6 code was recorded. |

Maternal-exposure codes P04.14, P04.16, and P04.41 are not used in these analyses. 

## Demographic variables

| Variable | Meaning |
|---|---|
| `age_years` | NCHS detail age converted to years; values outside 0–125 are missing. |
| `sex` | Male or Female after harmonizing numeric and character annual codes. |
| `race4` | White, Black, Asian/Pacific Islander (`API`), or Other. Other includes AIAN, multiracial, and unknown values. |
| `hisp3` | Non-Hispanic, Hispanic, or Unknown. Hispanic origin is separate from race. |

## Cause-of-death condition indicators

Each variable below equals 1 when at least one qualifying ICD-10 code appears as the underlying cause or in an available record-axis field and 0 otherwise. These indicators describe recorded causes of death, not a complete medical history.

| Variable | Condition domain |
|---|---|
| `alcohol_any` | F10, K70, X45, or T51 alcohol-attributable code. |
| `cardiac_hypertension` | Prespecified heart disease, selected vascular disease, or hypertension code. |
| `cancer` | Prespecified malignant-neoplasm or HIV-related neoplasm code. |
| `cerebrovascular` | I60–I69 cerebrovascular disease. |
| `diabetes_metabolic` | E10–E14 diabetes mellitus. |
| `hepatobiliary_pancreatic` | Prespecified liver, biliary, or pancreatic disease. |
| `renal` | Prespecified glomerular, acute/chronic kidney, or other renal disease. |
| `viral_hepatitis` | B15–B19 viral hepatitis. |
| `hiv_aids` | B20–B24 HIV/AIDS. |
| `substance_use_or_mental_health_disorder` | F10–F19 mental and behavioral disorders due to psychoactive substance use. |
| `accidental_injury_other_than_drug_overdose` | Prespecified external injury code, excluding drug-overdose intent codes; set to 0 for suicide or homicide manner of death. |

## Four modeled outcomes

| Variable | Construction |
|---|---|
| `cardiovascular` | Copy of `cardiac_hypertension`. |
| `cerebrovascular` | I60–I69 indicator described above. |
| `other_medical` | Equals 1 when any of cancer, diabetes/metabolic, hepatobiliary/pancreatic, renal, viral hepatitis, or HIV/AIDS equals 1. |
| `no_additional` | Equals 1 when cardiovascular, cerebrovascular, and other medical all equal 0. |

The exact ICD-10 regular expressions are documented directly above `condition_patterns` in `scripts/01_build_analysis_data.R`.
