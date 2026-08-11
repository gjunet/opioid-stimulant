# Derived Analysis Data Dictionary

Script 01 creates a local death-record-level RDS file with one row per death for all cases meeting the inclusive opioid or stimulant definition. The repository does not include mortality microdata. Researchers must obtain the annual National Center for Health Statistics Multiple Cause of Death files from the [CDC](https://www.cdc.gov/nchs/data_access/vitalstatsonline.htm) or another authorized source. The [NBER mortality-data documentation](https://www.nber.org/research/data/mortality-data-vital-statistics-nchs-multiple-cause-death-data) describes the annual source files and field coding.

## How ICD-10 categories are matched

The analysis removes punctuation and spaces from each ICD-10 value before matching. A category shown below without a decimal point includes every descendant subcode because the program matches the beginning of the normalized code. For example, `C50` includes `C50.0` through `C50.9`; `C00-C97` includes all malignant-neoplasm categories and their descendant subcodes. A single subcode such as `B21.2` includes only that subcode and any values recorded beneath it. The program searches the underlying cause and every available record-axis multiple-cause field (`recaxis1` through `recaxis20`).

The tables below reproduce the rules in `scripts/01_build_analysis_data.R`. They are the auditable specification of the implemented definitions, not a general list of every ICD-10 code that could describe a condition.

## Identification and exposure variables

| Variable | Meaning |
|---|---|
| `year` | Calendar year of death, 2015–2024. |
| `comparison_group` | Three mutually exclusive groups under the inclusive definition: `stimulant_no_opiate`, `opiate_stimulant`, or `opiate_no_stimulant`. |
| `restricted_group` | The same three groups reconstructed from poisoning and toxicologic blood-finding codes only. Missing means the death had an inclusive F11/F14/F15 code but no qualifying restricted code. |
| `opioid` | At least one F11, T40.0–T40.4, T40.6, or R78.1 code was recorded. |
| `cocaine` | At least one F14, T40.5, or R78.2 code was recorded. |
| `psychostimulant` | At least one F15 or T43.6 code was recorded. |

### Drug-involvement code hierarchy

| Drug class | Inclusive definition | Restricted definition | Notes |
|---|---|---|---|
| Opioid | `F11` (mental and behavioral disorders due to opioid use); `T40.0-T40.4` and `T40.6` (poisoning by, adverse effect of, or underdosing of the listed narcotics and related substances); `R78.1` (finding of opiate drug in blood). | `T40.0-T40.4`, `T40.6`, and `R78.1`. | `T40.5` is excluded from opioid involvement and assigned to stimulant involvement. |
| Cocaine | `F14` (mental and behavioral disorders due to cocaine use); `T40.5` (cocaine); `R78.2` (finding of cocaine in blood). | `T40.5` and `R78.2`. | Cocaine is one component of the broad stimulant definition. |
| Other psychostimulant | `F15` (mental and behavioral disorders due to other stimulant use, including caffeine); `T43.6` (psychostimulants with abuse potential). | `T43.6`. | This component is analyzed separately from cocaine in Sensitivity Analysis 2. |

The inclusive definition counts all descendant subcodes under `F11`, `F14`, and `F15`, not only selected fourth-character diagnoses. The restricted definition excludes every `F11`, `F14`, and `F15` code. Maternal-exposure codes `P04.14`, `P04.16`, and `P04.41` are excluded from both definitions.

## Demographic variables

| Variable | Meaning |
|---|---|
| `age_years` | NCHS detail age converted to years; values outside 0–125 are missing. |
| `sex` | Male or Female after harmonizing numeric and character annual codes. |
| `race4` | White, Black, Asian/Pacific Islander (`API`), or Other. Other includes AIAN, multiracial, and unknown values. |
| `hisp3` | Non-Hispanic, Hispanic, or Unknown. Hispanic origin is separate from race. |

## Cause-of-death condition indicators

Each variable below equals 1 when at least one qualifying ICD-10 code appears as the underlying cause or in an available record-axis field and 0 otherwise. These indicators describe conditions recorded on the death certificate. They do not represent a complete medical history, establish disease onset, or imply that drug involvement caused the condition.

| Condition indicator | ICD-10 categories implemented in the code | Category hierarchy represented |
|---|---|---|
| Alcohol involvement | `F10`; `K70`; `X45`; `T51` | Mental and behavioral disorders due to alcohol use; alcoholic liver disease; accidental poisoning by and exposure to alcohol; toxic effect of alcohol. All descendant subcodes are included. |
| Cardiac disease and hypertension | `I00-I09`; `I10-I13`; `I15`; `I20-I22`; `I24-I28`; `I30-I31`; `I33-I38`; `I40`; `I42-I51`; `I70-I71` | Rheumatic heart disease; selected hypertensive disease; ischemic heart disease; pulmonary heart disease and pulmonary-circulation disorders; selected pericardial, endocardial, valvular, myocardial, conduction, rhythm, heart-failure, and other heart-disease categories; atherosclerosis; aortic aneurysm and dissection. The definition does **not** include every code in Chapter IX. In particular, `I14`, `I16-I19`, `I23`, `I29`, `I32`, `I39`, `I41`, `I52-I59`, and `I72-I99` are not matched by this indicator. |
| Cancer | `B21.2`; `B21.8`; `C00-C97` | HIV disease resulting in other non-Hodgkin lymphoma (`B21.2`) or other malignant neoplasms (`B21.8`), plus the full malignant-neoplasm block: lip/oral cavity/pharynx (`C00-C14`); digestive organs (`C15-C26`); respiratory and intrathoracic organs (`C30-C39`); bone and articular cartilage (`C40-C41`); melanoma and other skin malignancies (`C43-C44`); mesothelial and soft tissue (`C45-C49`); breast (`C50`); female genital organs (`C51-C58`); male genital organs (`C60-C63`); urinary tract (`C64-C68`); eye, brain, and other central nervous system (`C69-C72`); thyroid and other endocrine glands (`C73-C75`); ill-defined, secondary, and unspecified sites (`C76-C80`); lymphoid, hematopoietic, and related tissue (`C81-C96`); and independent multiple-site malignancies (`C97`). |
| Cerebrovascular disease | `I60-I69` | The full cerebrovascular-disease block, including nontraumatic hemorrhage, cerebral infarction, other cerebrovascular diseases, and sequelae. |
| Diabetes | `E10-E14` | Type 1, type 2, malnutrition-related, other specified, and unspecified diabetes mellitus, including all descendant complication subcodes. No other endocrine, nutritional, or metabolic categories are included. |
| Hepatobiliary and pancreatic disease | `K70`; `K72-K74`; `K76`; `K80-K82`; `K85-K86` | Alcoholic liver disease; hepatic failure; chronic hepatitis; fibrosis/cirrhosis; other liver disease; cholelithiasis, cholecystitis, and other gallbladder disease; acute pancreatitis and other pancreatic disease. `K71`, `K75`, `K77-K79`, `K83-K84`, and `K87-K93` are not included. |
| Renal disease | `N00-N07`; `N17-N19`; `N25-N28` | Glomerular diseases; acute kidney failure and chronic kidney disease; unspecified kidney failure; disorders resulting from impaired renal tubular function; contracted kidney; small kidney of unknown cause; and other kidney/ureter disorders. `N08-N16`, `N20-N24`, and `N29-N99` are not included. |
| Viral hepatitis | `B15-B19` | The full viral-hepatitis block: hepatitis A, B, other acute viral hepatitis, chronic viral hepatitis, and unspecified viral hepatitis. |
| HIV/AIDS | `B20-B24` | The full human immunodeficiency virus disease block. |
| Substance-use or mental-health disorder | `F10-F19` | Mental and behavioral disorders due to psychoactive substance use: alcohol, opioids, cannabinoids, sedatives/hypnotics, cocaine, other stimulants, hallucinogens, tobacco, volatile solvents, and multiple/other substances. This descriptive indicator includes the same `F11`, `F14`, and `F15` families used in the inclusive drug definition. |
| Other accidental injury | `V00-V99`; `W00-W99`; `X00-X39`; `X60-X99`; `Y00-Y09`; `Y15-Y99` | Transport accidents; other external causes of accidental injury; exposure to smoke, fire, flames, heat, hot substances, forces of nature, and other specified accidental threats; intentional self-harm and assault categories; events of undetermined intent outside `Y10-Y14`; legal intervention/war; complications of medical and surgical care; and sequelae/supplementary external-cause factors. After code matching, the indicator is reset to 0 when the NCHS manner-of-death field identifies suicide or homicide. The implemented pattern excludes all `X40-X59` and `Y10-Y14`, not only drug-poisoning categories. |

The external-injury indicator is named for its analytic purpose, but the last row states the exact implemented rule. Readers should use those ranges, together with the manner-of-death override, when assessing whether the definition includes the external causes they expect.

### R variable names for the condition indicators

- Alcohol involvement: `alcohol_any`
- Cardiac disease and hypertension: `cardiac_hypertension`
- Cancer: `cancer`
- Cerebrovascular disease: `cerebrovascular`
- Diabetes: `diabetes_metabolic`
- Hepatobiliary and pancreatic disease: `hepatobiliary_pancreatic`
- Renal disease: `renal`
- Viral hepatitis: `viral_hepatitis`
- HIV/AIDS: `hiv_aids`
- Substance-use or mental-health disorder: `substance_use_or_mental_health_disorder`
- Other accidental injury: `accidental_injury_other_than_drug_overdose`

## Four modeled outcomes

| Variable | Construction |
|---|---|
| `cardiovascular` | Copy of `cardiac_hypertension`. |
| `cerebrovascular` | I60–I69 indicator described above. |
| `other_medical` | Equals 1 when any of cancer, diabetes/metabolic, hepatobiliary/pancreatic, renal, viral hepatitis, or HIV/AIDS equals 1. |
| `no_additional` | Equals 1 when cardiovascular, cerebrovascular, and other medical all equal 0. |

The exact ICD-10 regular expressions are documented directly above `condition_patterns` in `scripts/01_build_analysis_data.R`.
