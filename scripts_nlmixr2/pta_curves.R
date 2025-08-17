#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readr)
  library(ggplot2)
  library(glue)
})

expo <- readr::read_rds("results/exposure_summary.rds")

# ---- Choose ONE path below ----
use_ctr_trough <- FALSE

# Option A: Ctrough target
Ctrough_target <- 2.0  # mg/L (edit)
pta_A <- expo |>
  mutate(hit = Ctrough_tau >= Ctrough_target) |>
  group_by(dose_mg, ii_h) |>
  summarise(PTA = mean(hit), .groups = "drop") |>
  mutate(target = glue("Ctrough ≥ {Ctrough_target}"))

# Option B: AUC/MIC target (antibacterial)
MIC_values <- 2^(seq(-3, 0, by = 0.1))  # e.g., 0.125..8 mg/L
target_ratio <- 125  # e.g., AUC/MIC ≥ 125
pta_B <- tidyr::crossing(expo, MIC = MIC_values) |>
  mutate(hit = (AUC_tau / MIC) >= target_ratio) |>
  group_by(dose_mg, ii_h, MIC) |>
  summarise(PTA = mean(hit), .groups = "drop") |>
  mutate(target = glue("AUC/MIC ≥ {target_ratio}")) |> 
  filter(ii_h == 12)

if (use_ctr_trough) {
  plt <- ggplot(pta_A, aes(x = factor(dose_mg), y = PTA, group = factor(ii_h))) +
    geom_point() + geom_line(aes(linetype = factor(ii_h))) +
    scale_y_continuous(labels = scales::percent_format()) +
    labs(x = "Dose (mg)", y = "PTA", linetype = "Interval (h)",
         title = "PTA by Dose and Interval",
         subtitle = unique(pta_A$target)) +
    theme_bw()
  #out <- "results/pta_ctr_trough.png"
  #ggsave(out, plt, width = 7, height = 4.5, dpi = 300)
  # cat(glue("✅ Wrote PTA plot: {out}\n"))
} else {
  plt <- ggplot(pta_B, aes(x = MIC, y = PTA, color = factor(dose_mg))) +
    geom_line() + geom_point() +
    #facet_wrap(~ii_h, labeller = label_both) +
    #scale_x_log10() +
    scale_y_continuous(labels = scales::percent_format()) +
    labs(x = "MIC (mg/L)", y = "PTA", color = "Dose (mg)",
         title = "PTA vs MIC",
         subtitle = unique(pta_B$target)) +
    theme_bw()
  #out <- "results/pta_auc_mic.png"
  #ggsave(out, plt, width = 8, height = 5, dpi = 300)
  #cat(glue("✅ Wrote PTA plot: {out}\n"))
}
plt
