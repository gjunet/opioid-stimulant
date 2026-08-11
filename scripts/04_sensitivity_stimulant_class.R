# =============================================================================
# SCRIPT 04: SENSITIVITY 2, COCAINE AND OTHER PSYCHOSTIMULANTS SEPARATELY
#
# The primary analysis combines cocaine and other psychostimulants into one
# stimulant category. This post hoc analysis checks whether that pooling hides
# different patterns for the two stimulant classes.
# =============================================================================

# ---- 1. Load the inclusive cohort and initialize output containers -----------
source(file.path(publication_root,"R","helpers.R")); log_msg <- new_logger("04_sensitivity_stimulant_class.log")
d <- readRDS(analytic_path); results <- list(); samples <- list(); diagnostics <- list()
for(cls in c("Cocaine","Other psychostimulant")) {
  # ---- 2. Construct a nonoverlapping class-specific cohort -------------------
  # The cocaine analysis excludes every death with other-psychostimulant
  # involvement. The other-psychostimulant analysis excludes every death with
  # cocaine involvement. This prevents deaths involving both stimulant classes
  # from appearing in either class-specific reference group.
  if(cls=="Cocaine") { dd <- d[!d$psychostimulant,]; present <- dd$cocaine } else { dd <- d[!d$cocaine,]; present <- dd$psychostimulant }

  # Recreate the same three conceptual groups within the selected stimulant
  # class: class-specific stimulant-only, opioid plus that stimulant class, and
  # opioid-only. Records with neither the selected stimulant nor opioids are NA.
  raw_group <- ifelse(!dd$opioid & present,"stimulant_no_opiate",ifelse(dd$opioid & present,"opiate_stimulant",ifelse(dd$opioid & !present,"opiate_no_stimulant",NA_character_)))
  x <- prepare_model_data(dd,raw_group); counts <- table(x$grp)
  samples[[cls]] <- data.frame(stimulant_class=cls,total_n=nrow(x),stimulant_only_n=unname(counts["Stimulant-only"]),opioid_stimulant_n=unname(counts["Opioid-stimulant"]),opioid_only_n=unname(counts["Opioid-only"]))
  terms <- c("grpOpioid-stimulant","grpOpioid-only"); labels <- c(paste0("Opioid-",tolower(cls)," vs ",cls,"-only"),paste0("Opioid-only vs ",cls,"-only"))
  # ---- 3. Fit the common adjusted model for all four outcomes ----------------
  # The reference group changes with the analysis: cocaine-only for the cocaine
  # models and other-psychostimulant-only for the other-psychostimulant models.
  for(o in names(outcomes)) { fit <- fit_adjusted(o,x); z <- extract_terms(fit,terms,labels,o); z$stimulant_class <- cls; results[[paste(cls,o)]] <- z; diagnostics[[paste(cls,o)]] <- model_diagnostic(fit,o,paste("Stimulant class",cls)) }
}
# ---- 4. Write sample counts, model estimates, and diagnostics ----------------
write_table(do.call(rbind,samples),"sensitivity2_stimulant_class_samples.csv")
write_table(do.call(rbind,results),"sensitivity2_stimulant_class_models.csv")
write_diag(do.call(rbind,diagnostics),"sensitivity2_model_diagnostics.csv")
log_msg("Separate stimulant-class sensitivity complete")
