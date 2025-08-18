# placeholder_model_comp_with_formulas.R
# Fabricated model comparison using AIC/BIC formulas; last model is best.

# ---- Packages ----
library(tibble)
library(dplyr)
library(gt)

# ---- Inputs you can tweak ----
n_obs <- 480  # number of observations used for BIC (e.g., subjects × samples)

models <- tibble::tibble(
  Model = c(
    "1-compartment",
    "2-compartment",
    "3-compartment",
    "3C + 1 transit",
    "3C + 2 transit",
    "3C + 3 transit + IIVs (best)"
  ),
  k = c(3, 3, 3, 4, 4, 6)  # per your instructions
)

# ---- Fabricate large log-likelihoods (last model best) ----
# Make sure increments overcome the higher-parameter penalty for BIC too.
loglik_vals <- c(11952.5, 12220.7, 12437.1, 
                 12565.4, 12615.3, 12759.5)

tbl <- models |>
  mutate(
    `Log-likelihood` = loglik_vals,
    AIC = -2 * `Log-likelihood` + 2 * k,
    BIC = -2 * `Log-likelihood` + k * log(n_obs)
  ) |>
  mutate(
    `ΔAIC` = AIC - min(AIC),
    `ΔBIC` = BIC - min(BIC)
  ) |>
  select(Model, `ΔAIC`, `ΔBIC`, `Log-likelihood`) |>
  # lock in the display order
  mutate(Model = factor(Model, levels = Model)) |>
  arrange(Model)

# ---- Build gt table ----
gt_tbl <- tbl |>
  gt() |>
  tab_header(
    title = "Model Comparison",
  ) |>
  fmt_number(
    columns = c(`ΔAIC`, `ΔBIC`),
    decimals = 1
  ) |>
  fmt_number(
    columns = c(`Log-likelihood`),
    decimals = 0,
    use_seps = TRUE
  ) |>
  tab_style(
    style = list(
      cell_fill(color = "#F0FFF0"),
      cell_text(weight = "bold")
    ),
    locations = cells_body(rows = Model == "3C + 3 transit + IIVs (best)")
  ) |>
  cols_label(
    Model = "Model",
    `ΔAIC` = "ΔAIC",
    `ΔBIC` = "ΔBIC",
    `Log-likelihood` = "Log-likelihood"
  ) 

gt_tbl

# Optionally save to HTML:
# gtsave(gt_tbl, "placeholder_model_comparison_aic_bic.html")
