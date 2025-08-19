#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readr)
  library(ggplot2)
  library(glue)
  library(Hmisc)   # binconf()
})

# =============================================================================
# User settings
# =============================================================================
in_rds         <- "results/exposure_summary.rds"
use_ctr_trough <- FALSE              # TRUE = Option A; FALSE = Option B
ctr_trough_thr <- 0.1               # mg/L, if using Option A

# Option B (AUC/MIC target) settings:
MIC_values   <- 2^(seq(-3, 0, by = 0.1))  # 0.125 .. 1 mg/L (dense grid)
target_ratio <- 125                       # AUC/MIC threshold
keep_intervals <- NULL                      # set to e.g. c(12, 24) or NULL to keep all

# Outputs
out_png <- if (use_ctr_trough) "results/pta_ctr_trough_ci.png" else "results/pta_auc_mic_ci.png"
out_csv <- sub("\\.png$", ".csv", out_png)
out_rds <- sub("\\.png$", ".rds", out_png)

# =============================================================================
# Load data
# =============================================================================
expo <- readr::read_rds(in_rds)

# Graceful check
need_cols_A <- c("dose_mg", "ii_h", "Ctrough_tau")
need_cols_B <- c("dose_mg", "ii_h", "AUC_tau")
has_cols_A <- all(need_cols_A %in% names(expo))
has_cols_B <- all(need_cols_B %in% names(expo))
if (use_ctr_trough && !has_cols_A) stop(glue("Missing columns for Option A: {toString(setdiff(need_cols_A, names(expo)))}"))
if (!use_ctr_trough && !has_cols_B) stop(glue("Missing columns for Option B: {toString(setdiff(need_cols_B, names(expo)))}"))

# Keep selected intervals, if specified
if (!is.null(keep_intervals)) {
  expo <- expo |> filter(ii_h %in% keep_intervals)
}

#dir.create("results", showWarnings = FALSE)

# =============================================================================
# Helper to add Wilson CI columns
# =============================================================================
add_wilson_ci <- function(df, hits_col = "hits", n_col = "n") {
  ci <- Hmisc::binconf(df[[hits_col]], df[[n_col]], method = "wilson")
  df |>
    mutate(
      PTA      = .data[[hits_col]] / .data[[n_col]],
      PTA_lo   = pmax(ci[, "Lower"], 0),
      PTA_hi   = pmin(ci[, "Upper"], 1)
    )
}

# =============================================================================
# OPTION A: Ctrough target
# =============================================================================
if (use_ctr_trough) {
  pta_A <- expo |>
    mutate(hit = Ctrough_tau >= ctr_trough_thr) |>
    group_by(dose_mg, ii_h) |>
    summarise(hits = sum(hit), n = dplyr::n(), .groups = "drop") |>
    add_wilson_ci() |>
    mutate(
      interval = factor(ii_h, levels = sort(unique(ii_h)), labels = paste0(sort(unique(ii_h)), " h")),
      dose_f   = factor(dose_mg, levels = sort(unique(dose_mg))),
      target   = glue("Ctrough ≥ {ctr_trough_thr} mg/L")
    )

  plt <- ggplot(pta_A, aes(x = dose_f, y = PTA, group = interval, linetype = interval)) +
    geom_errorbar(aes(ymin = PTA_lo, ymax = PTA_hi), width = 0.15, alpha = 0.6) +
    geom_point(size = 2) +
    geom_line() +
    geom_hline(yintercept = 0.5, colour = "red") +
    scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
    labs(
      x = "Dose (mg)", y = "P(Target Attainment)",
      title = "PTA by Dose with 95% Binomial CI",
      subtitle = unique(pta_A$target),
      linetype = "Interval (h)"
    ) +
    theme_bw() +
    theme(legend.position = "bottom")

  # readr::write_rds(pta_A, out_rds)
  # readr::write_csv(pta_A, out_csv)
  # ggsave(out_png, plt, width = 7.5, height = 4.8, dpi = 300)
  # cat(glue("✅ Wrote: {out_png}\n"))

# =============================================================================
# OPTION B: AUC/MIC target (with ribbons/bands)
# =============================================================================
} else {
  # Build MIC cross-product and compute hits per (dose, interval, MIC)
  pta_B <- tidyr::crossing(expo, MIC = MIC_values) |>
    mutate(hit = (AUC_tau / MIC) >= target_ratio) |>
    group_by(dose_mg, ii_h, MIC) |>
    summarise(hits = sum(hit), n = dplyr::n(), .groups = "drop") |>
    add_wilson_ci() |>
    mutate(
      dose_f   = factor(dose_mg, levels = sort(unique(dose_mg))),
      interval = factor(ii_h, levels = sort(unique(ii_h))),
      target   = glue("AUC/MIC ≥ {target_ratio}")
    ) |>
    arrange(ii_h, dose_mg, MIC)

  # Single-line view (like your current plot):
  # - one panel (since we filtered to ii_h == 12 above)
  # - colored by dose
  plt <- ggplot(pta_B, aes(x = MIC, y = PTA, color = dose_f)) +
    geom_ribbon(aes(ymin = PTA_lo, ymax = PTA_hi, group = dose_f, fill = dose_f),
                alpha = 0.18, color = NA) +
    geom_line(size = 0.9) +
    geom_point(size = 1.4) +
    geom_hline(yintercept = 0.5, linetype = "dashed") +
    geom_vline(xintercept = 0.1, linetype = "dashed") +
    # If you want log-x, uncomment:
    # scale_x_log10() +
    scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
    labs(
      x = "MIC (mg/L)", y = "P(Target Attainment)",
      color = "Dose (mg)", fill = "Dose (mg)",
      title = "PTA vs MIC with 95% Binomial CI (Wilson)",
      subtitle = unique(pta_B$target)
    ) +
    theme_bw() +
    theme(legend.position = "bottom")

  #readr::write_rds(pta_B, out_rds)
  #readr::write_csv(pta_B, out_csv)
  #ggsave(out_png, plt, width = 9, height = 5.2, dpi = 300)
  #cat(glue("✅ Wrote: {out_png}\n"))
}

# Also print to device if interactive context (keeps your current behavior)
plt
