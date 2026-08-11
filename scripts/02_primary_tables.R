# =============================================================================
# SCRIPT 02: PRIMARY MANUSCRIPT TABLES
#
# This script describes the inclusive 2015-2024 cohort and estimates the main
# associations between drug-involvement group and each medical-condition
# outcome. It creates manuscript Table 1, Appendix Tables S1-S7, an ICD-10 drug
# code table, and model diagnostics.
# =============================================================================

# ---- 1. Load the derived cohort and apply common factor/reference coding -----
source(file.path(publication_root,"R","helpers.R")); log_msg <- new_logger("02_primary_tables.log")
d <- readRDS(analytic_path); x <- prepare_model_data(d,d$comparison_group)

# ---- 2. Create Table 1 and Appendix Table S1: sample characteristics ---------
# For binary and categorical characteristics, cells show n (column %). Age is
# shown as mean (SD). Pearson chi-square tests compare groups for categorical
# variables; one-way ANOVA compares mean age. These P values are descriptive
# omnibus comparisons and are not the adjusted regression results.
groups <- levels(x$grp); labels <- c(cardiovascular="Cardiovascular",cerebrovascular="Cerebrovascular",other_medical="Other medical",no_additional="No additional")
np <- function(v,idx) sprintf("%s (%.1f%%)",format(sum(v[idx]),big.mark=","),100*mean(v[idx]))
rows <- list(); add <- function(characteristic,values,p="") rows[[length(rows)+1]] <<- data.frame(Characteristic=characteristic,`Stimulant-only`=values[1],`Opioid-stimulant`=values[2],`Opioid-only`=values[3],p_value=p,check.names=FALSE)
idx <- lapply(groups,function(g)x$grp==g)
add("N",vapply(idx,function(i)format(sum(i),big.mark=","),character(1)))
for(o in names(labels)) add(labels[o],vapply(idx,function(i)np(x[[o]],i),character(1)),fmt_p(chisq.test(table(x$grp,x[[o]]))$p.value))
add("Age, mean (SD)",vapply(idx,function(i)sprintf("%.1f (%.1f)",mean(x$age_years[i],na.rm=TRUE),sd(x$age_years[i],na.rm=TRUE)),character(1)),fmt_p(summary(aov(age_years~grp,data=x))[[1]][["Pr(>F)"]][1]))
for(spec in list(c("sex","Sex"),c("race4","Race"),c("hisp3","Hispanic origin"))) {
  variable <- spec[1]; heading <- spec[2]; p <- fmt_p(chisq.test(table(x$grp,x[[variable]]))$p.value)
  for(i in seq_along(levels(x[[variable]]))) {
    level <- levels(x[[variable]])[i]; indicator <- as.integer(x[[variable]]==level)
    add(paste0(heading,": ",level),vapply(idx,function(j)np(indicator,j),character(1)),if(i==1)p else "")
  }
}
table1 <- do.call(rbind,rows); write_table(table1,"table1_chang_2015_2024.csv"); write_table(table1,"S1_descriptive_2015_2024.csv")

# ---- 3. Fit unadjusted and adjusted logistic regression models ---------------
# One model is fit for each binary outcome. The unadjusted models contain drug
# group only. The adjusted models add age (per 10 years), sex, race, and
# Hispanic origin. Stimulant-only is the reference, so the two exposure terms
# compare opioid-stimulant and opioid-only deaths with stimulant-only deaths.
unadjusted <- list(); adjusted <- list(); diagnostics <- list(); pair_files <- c(cardiovascular="S4_pairwise_cardiovascular.csv",cerebrovascular="S5_pairwise_cerebrovascular.csv",other_medical="S6_pairwise_other_medical.csv",no_additional="S7_pairwise_no_additional.csv")
terms <- c("grpOpioid-stimulant","grpOpioid-only"); contrasts <- c("Opioid-stimulant vs Stimulant-only","Opioid-only vs Stimulant-only")
for(o in names(outcomes)) {
  fu <- glm(as.formula(paste(o,"~ grp")),data=x,family=binomial()); fa <- fit_adjusted(o,x)
  unadjusted[[o]] <- extract_terms(fu,terms,contrasts,o); adjusted[[o]] <- extract_terms(fa,terms,contrasts,o)
  diagnostics[[o]] <- model_diagnostic(fa,o,"Primary adjusted")
  # Use the fitted coefficient covariance matrix to obtain all three pairwise
  # drug-group contrasts. The third contrast subtracts the opioid-stimulant
  # coefficient from the opioid-only coefficient. Bonferroni correction
  # multiplies each two-sided P value by three because three pairs are tested.
  b <- coef(fa); V <- vcov(fa); est <- c(b[terms[1]],b[terms[2]],b[terms[2]]-b[terms[1]])
  se <- c(sqrt(V[terms[1],terms[1]]),sqrt(V[terms[2],terms[2]]),sqrt(V[terms[1],terms[1]]+V[terms[2],terms[2]]-2*V[terms[1],terms[2]]))
  names(est) <- c(contrasts,"Opioid-only vs Opioid-stimulant")
  p <- pmin(1,2*pnorm(-abs(est/se))*3)
  pair <- data.frame(contrast=names(est),OR=exp(est),CI_low=exp(est-1.96*se),CI_high=exp(est+1.96*se),p_bonferroni=p)
  write_table(pair,pair_files[o])
}
format_model <- function(z) { z <- do.call(rbind,z); data.frame(outcome=z$outcome,comparison=z$contrast,OR=z$odds_ratio,CI_low=z$or_ci_low,CI_high=z$or_ci_high,p_value=z$p_value) }
write_table(format_model(unadjusted),"S2_unadjusted_2015_2024.csv")
write_table(format_model(adjusted),"S3_adjusted_2015_2024.csv")
write_diag(do.call(rbind,diagnostics),"primary_model_diagnostics.csv")

# ---- 4. Document drug-involvement codes used by both definitions -------------
# "Inclusive only" identifies F-code families removed from the restricted
# sensitivity definition. "Inclusive and restricted" identifies poisoning and
# toxicologic blood-finding codes shared by both definitions. The three P04
# codes are displayed to document their deliberate exclusion from both.
codes <- data.frame(domain="Drug involvement",drug_class=c(rep("Opioid",8),rep("Stimulant",5),rep("Excluded from both",3)),code=c("F11 (all subcodes)","T40.0","T40.1","T40.2","T40.3","T40.4","T40.6","R78.1","F14 (all subcodes)","F15 (all subcodes)","T40.5","T43.6","R78.2","P04.14","P04.16","P04.41"),definition_status=c("Inclusive only",rep("Inclusive and restricted",7),"Inclusive only","Inclusive only",rep("Inclusive and restricted",3),rep("Excluded from both",3)))
write_table(codes,"eTable_icd10_drug_codes.csv")
log_msg("Primary tables complete; N=",format(nrow(d),big.mark=","))
