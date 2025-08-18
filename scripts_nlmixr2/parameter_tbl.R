
fx <- tidy(fit, effects = "fixed")
rp <- tidy(fit, effects = "ran_pars") |> dplyr::filter(group=="ID")  
sg <- tidy(fit, effects = "ran_pars") |> dplyr::filter(group!="ID")  

# # --- Fixed effects: back-transform log-params ---
bt_map <- c(tcl = "CL", tvc = "Vc", tktr = "ktr")
fixed_tbl <- fx %>%
  transmute(
    term,
    estimate,
    se = dplyr::coalesce(std.error, NA_real_)
  ) %>%
  mutate(
    Parameter   = dplyr::recode(term, !!!bt_map, .default = term),
    Estimate_bt = if_else(term %in% names(bt_map), exp(estimate), estimate),
    SE_bt       = if_else(term %in% names(bt_map) & !is.na(se), exp(estimate) * se, se),
    LCI         = if_else(!is.na(SE_bt), Estimate_bt * exp(-1.96 * dplyr::coalesce(se, 0)), NA_real_),
    UCI         = if_else(!is.na(SE_bt), Estimate_bt * exp( 1.96 * dplyr::coalesce(se, 0)), NA_real_)
  ) %>%
  transmute(
    Component = "Fixed effects",
    Parameter,
    Estimate  = Estimate_bt,
    SE        = SE_bt,
    LCI, UCI,
    Notes     = NA_character_
  )

# --- Random effects: parse SDs (diag) and correlations (off-diag) ---
# broom.mixed typically returns columns like: effect, group, term, estimate, component (\"sd\", \"cor\", or \"sdcor\")
iiv_var <- NULL
iiv_cor <- NULL


# SDs -> convert to variances; add %CV note (SD ~ omega on ETA-scale)
iiv_var <- rp |> 
    transmute(
        Component = "Random effects (Ω)",
        Parameter = dplyr::case_when(
            term == "sd__eta.cl" ~ "IIV CL",
            term == "sd__eta.vc" ~ "IIV Vc",
            TRUE ~ term),
        Estimate  = estimate^2,  # variance
        SE        = NA_real_, LCI = NA_real_, UCI = NA_real_,
        Notes     = paste0("CV%≈ ", sprintf("%.1f", estimate * 100))
        )

# --- Residual (Σ) ---
sigma_tbl <- sg |> 
    transmute(
      Component = "Residual error (Σ)",
      Parameter = "Proportional",
      Estimate  = estimate,
      SE = NA_real_,
      Notes = NA_character_
    )

# --- Combine and render ---
all_params <- bind_rows(fixed_tbl, iiv_var, iiv_cor, sigma_tbl) %>%
  mutate(Component = factor(Component, levels = c("Fixed effects","Random effects (Ω)","Residual error (Σ)"))) %>%
  arrange(Component, desc(str_detect(Parameter, "^IIV\\(")), Parameter)