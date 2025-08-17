#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggdist)
  library(ggplot2)
  library(glue)
  library(stringr)
})

in_expo <- "results/exposure_summary.rds"
out_png <- "results/exposure_facets.png"

expo <- readr::read_rds(in_expo)

# Long format: choose the exposure metrics to plot
metrics <- c("AUC_tau", "Cmax_tau", "Ctrough_tau", "Cavg_tau")

expo_long <- expo |>
  select(id, dose_mg, ii_h, all_of(metrics)) |>
  pivot_longer(cols = all_of(metrics), names_to = "metric", values_to = "value") |>
  mutate(
    metric = factor(metric,
                    levels = c("AUC_tau", "Cmax_tau", "Ctrough_tau", "Cavg_tau"),
                    labels = c("AUC₀–τ", "Cmax (τ)", "Ctrough (τ)", "Cavg (τ)")),
    interval = factor(ii_h, levels = sort(unique(ii_h)), labels = paste0(sort(unique(ii_h)), " h")),
    dose_f = factor(dose_mg, levels = sort(unique(dose_mg)))
  )

# Plot: boxplots (per-subject distribution) + median line across doses


p<- (
  ggplot(expo_long |> mutate(obs = if_else(dose_f == 600, "yes", "no")), 
    aes(x = dose_f, y = value, colour = obs)) +
    geom_boxplot(outlier.size = 0.4, width = 0.7) +
    facet_grid(rows = vars(metric), cols = vars(interval), scales = "free_y") +
    scale_y_continuous(transform = scales::pseudo_log_trans()) +
    scale_colour_manual(values = c("yes" = "grey50", "no" = "black")) +
    labs(
      x = "Dose (mg)",
      y = "Value",
      title = "Exposure Metrics by Dose",
      subtitle = "Faceted by metric (rows) and dosing interval (columns)"
    ) +
    theme_bw() +
    theme(
      legend.position = "None",
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey90"),
      axis.title.y = element_text(margin = margin(r = 8)),
      axis.title.x = element_text(margin = margin(t = 8))
  )
)

dir.create("results", showWarnings = FALSE)
ggsave(out_png, p, width = 10, height = 8, dpi = 300)
cat(glue("✅ Wrote faceted exposure plot: {out_png}\n"))
