# =============================================================================
# SCRIPT 08: FIGURE 1 FOREST PLOTS
#
# This script recreates Figure 1 from the adjusted logistic regression models.
# It makes one version on the odds-ratio scale and one on the unexponentiated
# logit-coefficient scale. Point labels report the estimate and its 95%
# confidence interval. The figure has no title; the reference group, adjustment
# set, and data-access statement appear below the graph.
#
# The manuscript models use stimulant-only deaths as the reference group. For
# transparency, the script also creates companion figures with opioid-only
# deaths as the reference. Changing the reference group changes the contrasts,
# but it does not change model fit or the underlying fitted probabilities.
# =============================================================================

# ---- 1. Load the cohort and plotting software -------------------------------
source(file.path(publication_root, "R", "helpers.R"))
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Figure 1 requires ggplot2. Install it with install.packages('ggplot2').")
}
library(ggplot2)
log_msg <- new_logger("08_figure1_forest_plots.log")
d <- readRDS(analytic_path)

# ---- 2. Extract the two drug-group contrasts for each outcome ---------------
# Coefficients are taken directly from the same adjusted logistic models used
# for the manuscript tables. Wald 95% confidence intervals are shown on both
# scales; odds ratios and their intervals are exponentiated coefficients.
figure_estimates <- function(reference_group) {
  x <- prepare_model_data(d, d$comparison_group)
  x$grp <- relevel(x$grp, ref = reference_group)
  comparison_groups <- setdiff(levels(x$grp), reference_group)
  terms <- paste0("grp", comparison_groups)
  contrast_labels <- paste0(comparison_groups, " vs ", reference_group)

  rows <- lapply(names(outcomes), function(outcome) {
    fit <- fit_adjusted(outcome, x)
    extract_terms(fit, terms, contrast_labels, outcome)
  })
  estimates <- do.call(rbind, rows)
  estimates$outcome <- factor(estimates$outcome, levels = unname(outcomes))
  estimates$contrast <- factor(
    estimates$contrast,
    levels = rev(contrast_labels)
  )
  estimates
}

# ---- 3. Apply a common publication style -----------------------------------
# Color and point shape both identify the comparison so the figure remains
# interpretable in grayscale. Labels are aligned at the right edge of each
# panel, leaving the point and confidence interval unobstructed.
comparison_colors <- c("#0072B2", "#D55E00")
comparison_shapes <- c(16, 17)
data_caption <- paste0(
  "Point labels show the estimate (95% CI).\n",
  "Models adjusted for age, sex, race, and Hispanic origin.\n",
  "Data accessed through the National Center for Health Statistics Research Data Center."
)

base_figure <- function(estimates, reference_group) {
  ggplot(estimates, aes(y = contrast, color = contrast, shape = contrast)) +
    facet_wrap(~ outcome, ncol = 1) +
    scale_color_manual(values = comparison_colors, guide = "none") +
    scale_shape_manual(values = comparison_shapes, guide = "none") +
    labs(
      y = NULL,
      caption = paste0("Reference group: ", reference_group, " deaths.\n", data_caption)
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold"),
      plot.caption = element_text(size = 8.5, hjust = 0),
      plot.caption.position = "plot",
      plot.margin = margin(8, 12, 8, 8)
    )
}

# ---- 4. Draw the adjusted odds-ratio version --------------------------------
# The odds-ratio axis is logarithmic, with 1 marking no difference from the
# reference group. The point label reports OR (95% CI).
make_or_figure <- function(estimates, reference_group) {
  lower_limit <- min(estimates$or_ci_low) * 0.75
  upper_limit <- max(estimates$or_ci_high) * 2.7
  estimates$value_label <- sprintf(
    "%.2f (%.2f-%.2f)",
    estimates$odds_ratio, estimates$or_ci_low, estimates$or_ci_high
  )
  estimates$label_x <- ifelse(
    estimates$or_ci_high < 1,
    1.25,
    estimates$or_ci_high * 1.16
  )

  base_figure(estimates, reference_group) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
    geom_errorbar(
      aes(xmin = or_ci_low, xmax = or_ci_high),
      width = 0.18, linewidth = 0.7, orientation = "y"
    ) +
    geom_point(aes(x = odds_ratio), size = 2.6) +
    geom_text(
      aes(x = label_x, label = value_label),
      hjust = 0, color = "black", size = 3.1, show.legend = FALSE
    ) +
    scale_x_log10(limits = c(lower_limit, upper_limit)) +
    labs(x = "Adjusted odds ratio (log scale)")
}

# ---- 5. Draw the adjusted logit-coefficient version -------------------------
# These are the unexponentiated coefficients from the adjusted logistic model.
# Zero marks no difference from the reference group. The point label reports
# beta (95% CI) on the log-odds scale.
make_logit_figure <- function(estimates, reference_group) {
  plotted_min <- min(estimates$beta_ci_low)
  plotted_max <- max(estimates$beta_ci_high)
  span <- plotted_max - plotted_min
  lower_limit <- plotted_min - 0.08 * span
  upper_limit <- plotted_max + 0.58 * span
  estimates$value_label <- sprintf(
    "%.2f (%.2f to %.2f)",
    estimates$beta, estimates$beta_ci_low, estimates$beta_ci_high
  )
  estimates$label_x <- ifelse(
    estimates$beta_ci_high < 0,
    0.30,
    estimates$beta_ci_high + 0.08 * span
  )

  base_figure(estimates, reference_group) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_errorbar(
      aes(xmin = beta_ci_low, xmax = beta_ci_high),
      width = 0.18, linewidth = 0.7, orientation = "y"
    ) +
    geom_point(aes(x = beta), size = 2.6) +
    geom_text(
      aes(x = label_x, label = value_label),
      hjust = 0, color = "black", size = 3.1, show.legend = FALSE
    ) +
    scale_x_continuous(limits = c(lower_limit, upper_limit)) +
    labs(x = "Adjusted logit coefficient")
}

# ---- 6. Save primary and alternate-reference figures ------------------------
# PNG files are suitable for manuscript review; SVG files remain editable and
# scale without loss of resolution. A CSV records every plotted value.
save_figure_pair <- function(reference_group, suffix) {
  estimates <- figure_estimates(reference_group)
  or_plot <- make_or_figure(estimates, reference_group)
  logit_plot <- make_logit_figure(estimates, reference_group)

  or_stem <- paste0("Figure1_adjusted_odds_ratios", suffix)
  logit_stem <- paste0("Figure1_adjusted_logit_coefficients", suffix)
  ggsave(file.path(figure_dir, paste0(or_stem, ".png")), or_plot,
         width = 9, height = 7.5, dpi = 300)
  ggsave(file.path(figure_dir, paste0(or_stem, ".svg")), or_plot,
         width = 9, height = 7.5)
  ggsave(file.path(figure_dir, paste0(logit_stem, ".png")), logit_plot,
         width = 9, height = 7.5, dpi = 300)
  ggsave(file.path(figure_dir, paste0(logit_stem, ".svg")), logit_plot,
         width = 9, height = 7.5)

  plotted_values <- estimates[, c(
    "outcome", "contrast", "beta", "beta_ci_low", "beta_ci_high",
    "odds_ratio", "or_ci_low", "or_ci_high", "p_value", "model_n"
  )]
  write.csv(
    plotted_values,
    file.path(table_dir, paste0("Figure1_plotted_estimates", suffix, ".csv")),
    row.names = FALSE
  )
  log_msg("Saved OR and logit Figure 1 files for reference group: ", reference_group)
}

save_figure_pair("Stimulant-only", "")
save_figure_pair("Opioid-only", "_opioid_reference")
