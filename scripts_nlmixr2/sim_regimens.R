#!/usr/bin/env Rscript
rm(list = ls())

suppressPackageStartupMessages({
  library(rxode2)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readr)
  library(gt)
  library(glue)
  library(tibble)
  library(rlang)
  library(pracma)   # trapz
})

# -------------------------------
# MODEL (same structure as your script)
# -------------------------------
model <- rxode2({
  d/dt(A0)  = -ktr * A0;
  d/dt(TR1) =  ktr * A0 - ktr * TR1;
  d/dt(TR2) =  ktr * TR1 - ktr * TR2;
  d/dt(TR3) =  ktr * TR2 - ktr * TR3;

  d/dt(A1) = ktr * TR3 - (CL / Vc) * A1
             - Q2 / Vc * A1 + Q2 / Vp1 * A2
             - Q3 / Vc * A1 + Q3 / Vp2 * A3;
  d/dt(A2) = Q2 / Vc * A1 - Q2 / Vp1 * A2;
  d/dt(A3) = Q3 / Vc * A1 - Q3 / Vp2 * A3;

  Cp = A1 / Vc;
})

# -------------------------------
# Population + IIV generator
# -------------------------------
mk_subject_params <- function(n_subjects = 50, seed = 101) {
  set.seed(seed)

  theta <- c(
    CL = 15, Vc = 60, Vp1 = 30, Vp2 = 40, Q2 = 6.0, Q3 = 5.0, MTT = 10, NTR = 3
  )
  theta_wt <- 0.75
  ktr <- (theta["NTR"] + 1) / theta["MTT"]

  # weights ~ lognormal
  mean_wt <- 70
  cv_wt <- 0.3
  omega_wt <- sqrt(log(1 + cv_wt^2))
  log_mean <- log(mean_wt) - 0.5 * omega_wt^2
  WT <- rlnorm(n_subjects, meanlog = log_mean, sdlog = omega_wt)

  # IIV (lognormal on CL, Vc)
  cv_iiv <- c(CL = 0.3, Vc = 0.2)
  Omega <- diag(log(1 + cv_iiv^2))
  etas <- MASS::mvrnorm(n_subjects, mu = rep(0, 2), Sigma = Omega)

  tibble(
    id = seq_len(n_subjects),
    WT = WT,
    eta.cl = etas[, 1],
    eta.vc = etas[, 2]
  ) |>
    mutate(
      CL  = 15 * (WT / 70)^theta_wt * exp(eta.cl),
      Vc  = 60 * exp(eta.vc),
      Vp1 = 30, Vp2 = 40, Q2 = 6.0, Q3 = 5.0,
      MTT = 10, NTR = 3, ktr = ktr
    )
}

# -------------------------------
# Regimen grid
# - Edit doses/ii/duration_h to explore quickly
# -------------------------------
regimen_grid <- tidyr::crossing(
  dose_mg    = c(200, 400, 600, 800),
  ii_h       = c(12, 24, 48, 72),
  duration_h = c(72*2)   # simulate 6 days by default; adjust as needed
) |>
  mutate(taus = ii_h)

# -------------------------------
# Build events for a regimen
# -------------------------------
mk_events <- function(n_subjects, dose_mg, ii_h, duration_h) {
  n_doses <- floor(duration_h / ii_h)
  obs_times <- sort(unique(c(seq(0, duration_h, by = 6), seq(ii_h, duration_h, by = ii_h) - 1e-6)))

  doses <- tibble(
    id = 1:n_subjects,
    time = 0,
    amt = dose_mg,
    ii = ii_h,
    addl = max(n_doses - 1, 0L),
    cmt = 1,
    evid = 1
  )

  obs <- tibble(
    id = rep(1:n_subjects, each = length(obs_times)),
    time = rep(obs_times, times = n_subjects),
    amt = 0, ii = 0, addl = 0, cmt = 1, evid = 0
  )

  bind_rows(doses, obs) |> arrange(id, time, desc(evid))
}

# -------------------------------
# Exposure metrics per subject over [0, duration] and over last dosing interval [t_end-τ, t_end]
# -------------------------------
compute_exposures <- function(sim_df, tau_h, duration_h) {
  end <- duration_h
  start_ss <- max(end - tau_h, 0)

  # whole window
  by_id <- sim_df |>
    group_by(id) |>
    summarise(
      AUC_0_T = trapz(time, Cp),
      Cmax_0_T = max(Cp),
      Tmax_0_T = time[which.max(Cp)[1]],
      .groups = "drop"
    )

  # last interval ~ pseudo steady-state interval if duration is long enough
  win <- sim_df |> filter(time >= start_ss, time <= end)
  by_id_ss <- win |>
    group_by(id) |>
    summarise(
      AUC_tau = trapz(time, Cp),
      Cmax_tau = max(Cp),
      Ctrough_tau = Cp[which.max(time)],  # last timepoint in window
      Cavg_tau = AUC_tau / tau_h,
      .groups = "drop"
    )

  left_join(by_id, by_id_ss, by = "id")
}

# -------------------------------
# Main
# -------------------------------
n_subjects <- 50
subject_params <- mk_subject_params(n_subjects = n_subjects, seed = 101)

# run all regimens
results_list <- regimen_grid |>
  mutate(regimen_id = row_number()) |>
  pmap(function(dose_mg, ii_h, duration_h, taus, regimen_id) {
    ev <- mk_events(n_subjects, dose_mg, ii_h, duration_h)
    sim <- rxSolve(
      model,
      params = subject_params |>
        dplyr::select(id, CL, Vc, Vp1, Vp2, Q2, Q3, ktr),
      events = ev,
      addDosing = TRUE
    ) |>
      as_tibble()

    expo <- compute_exposures(sim_df = sim, tau_h = ii_h, duration_h = duration_h) |>
      mutate(regimen_id = regimen_id, dose_mg = dose_mg, ii_h = ii_h, duration_h = duration_h)

    list(sim = sim |> mutate(regimen_id = regimen_id),
         expo = expo)
  })

# bind exposures across regimens
exposure_all <- map_dfr(results_list, ~.x$expo)

# (Optional) also keep a slim simulation sample for plotting
sim_sample <- map_dfr(results_list, ~.x$sim |> semi_join(
  tibble(id = sample(unique(.x$sim$id), min(8, length(unique(.x$sim$id))))), by = "id"
))

# save
out_expo <- "./results/exposure_summary.rds"
out_sim  <- "./results/sim_sample.rds"

write_rds(exposure_all, out_expo)
write_rds(sim_sample, out_sim)

cat(glue("✅ Saved exposures: {out_expo}\n"))
cat(glue("✅ Saved slim sims: {out_sim}\n"))
