#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(gt)
  library(glue)
})

in_expo <- "results/exposure_summary.rds"
exposure_all <- readr::read_rds(in_expo)

# Summary by regimen (median [IQR])
summ_tab <- exposure_all |>
  group_by(dose_mg, ii_h) |>
  summarise(
    n = dplyr::n(),
    AUC_tau_med = median(AUC_tau, na.rm = TRUE),
    AUC_tau_iqr = IQR(AUC_tau, na.rm = TRUE),
    Cmax_tau_med = median(Cmax_tau, na.rm = TRUE),
    Cmax_tau_iqr = IQR(Cmax_tau, na.rm = TRUE),
    Ctrough_tau_med = median(Ctrough_tau, na.rm = TRUE),
    Ctrough_tau_iqr = IQR(Ctrough_tau, na.rm = TRUE),
    Cavg_tau_med = median(Cavg_tau, na.rm = TRUE),
    Cavg_tau_iqr = IQR(Cavg_tau, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    AUC_tau = glue("{signif(AUC_tau_med,3)} [{signif(AUC_tau_med - 0.5*AUC_tau_iqr,3)}, {signif(AUC_tau_med + 0.5*AUC_tau_iqr,3)}]"),
    Cmax_tau = glue("{signif(Cmax_tau_med,3)} [{signif(Cmax_tau_med - 0.5*Cmax_tau_iqr,3)}, {signif(Cmax_tau_med + 0.5*Cmax_tau_iqr,3)}]"),
    Ctrough_tau = glue("{signif(Ctrough_tau_med,3)} [{signif(Ctrough_tau_med - 0.5*Ctrough_tau_iqr,3)}, {signif(Ctrough_tau_med + 0.5*Ctrough_tau_iqr,3)}]"),
    Cavg_tau = glue("{signif(Cavg_tau_med,3)} [{signif(Cavg_tau_med - 0.5*Cavg_tau_iqr,3)}, {signif(Cavg_tau_med + 0.5*Cavg_tau_iqr,3)}]")
  ) |>
  select(dose_mg, ii_h, n, AUC_tau, Cmax_tau, Ctrough_tau, Cavg_tau) |>
  arrange(ii_h, dose_mg)

tbl <- summ_tab |>
  gt() |>
  tab_header(title = "Exposure by Regimen",
             subtitle = "Median [±0.5·IQR] per subject across last dosing interval") |>
  cols_label(
    dose_mg = "Dose (mg)",
    ii_h = "Interval (h)",
    n = "N",
    AUC_tau = "AUC₀–τ",
    Cmax_tau = "Cmax",
    Ctrough_tau = "Ctrough",
    Cavg_tau = "Cavg"
  )

out_tbl <- "results/exposure_table.html"
gtsave(tbl, out_tbl)
cat(glue("✅ Wrote exposure table: {out_tbl}\n"))
