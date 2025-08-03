# requires ---
# tidyverse
# mrgsolve

# simulate individuals
simulate_individuals <- function(n_per_arm = 10, seed = 123) {
  set.seed(seed)

  regimens <- c("0.5_q2w", "2_q2w", "5_q2w", "10_q2w", "2_q3w", "5_q3w", "200_q3w")
  n_total <- length(regimens) * n_per_arm

  tibble(
    id = 1:n_total,
    regimen = rep(regimens, each = n_per_arm),
    sex = rbinom(n_total, 1, prob = 0.74),  # 1 = Male, 0 = Female
    wt = rlnorm(n_total, meanlog = log(65), sdlog = 0.25),
    age = rlnorm(n_total, meanlog = log(60), sdlog = 0.20)
  ) |>
    mutate(
      sex = if_else(sex == 1, "M", "F"),
      wt = round(wt, 1),
      age = round(age, 1)
    )
}

# simulate indvidual schedule
create_event_schedule <- function(individuals, sim_days = 180) {
  regimen_info <- tibble(
    regimen = c("0.5_q2w", "2_q2w", "5_q2w", "10_q2w", "2_q3w", "5_q3w", "200_q3w"),
    dose_mgkg = c(0.5, 2, 5, 10, 2, 5, NA),
    freq = c(14, 14, 14, 14, 21, 21, 21),
    is_flat = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE),
    dose_mg = c(NA, NA, NA, NA, NA, NA, 200)
  )

  individuals <- individuals |>
    left_join(regimen_info, by = "regimen") |>
    mutate(
      dose_amt = if_else(
        is_flat,
        dose_mg,
        dose_mgkg * wt
      ),
      n_doses = floor(sim_days / freq)
    )

  # Build event table for each individual
  all_events <- individuals |>
    rowwise() |>
    mutate(
      ev = list(
        ev(
          ID = id,
          amt = dose_amt,
          cmt = 1,
          ii = freq,
          addl = n_doses - 1,
          time = 0
        )
      )
    ) |>
    pull(ev) |> 
    map(as_tibble) |> 
    bind_rows()
  
  return(all_events)
}

# simulate individual parameter data 
generate_idata <- function(individuals, seed = 202) {
  set.seed(seed)

  N <- nrow(individuals)

  eta <- tibble(
    # PK parameters
    ETA_CL = rnorm(N, 0, sqrt(log(1 + (0.263)^2))),
    ETA_VC = rnorm(N, 0, sqrt(log(1 + (0.167)^2))),
    ETA_V2 = rnorm(N, 0, sqrt(log(1 + (0.747)^2))),
    ETA_V3 = rnorm(N, 0, sqrt(log(1 + (0.999)^2))),
    # PD parameters
    ETA_KIN = rnorm(N, 0, sqrt(log(1 + 0.3^2))),
    ETA_KOUT = rnorm(N, 0, sqrt(log(1 + 0.3^2))),
    ETA_EC50 = rnorm(N, 0, sqrt(log(1 + 0.4^2)))
  )

  individuals |>
    bind_cols(eta) |>
    mutate(
      # PK params as before
      CL = 0.153 * (wt / 65)^0.565 * exp(ETA_CL),
      VC = 3.05 * (wt / 65)^0.397 * exp(ETA_VC),
      Q2 = 0.74,
      V2 = 1.27 * exp(ETA_V2),
      Q3 = 0.092,
      V3 = 2.10 * exp(ETA_V3),
      # PD params
      KIN = 10 * exp(ETA_KIN),
      KOUT = 1 * exp(ETA_KOUT),
      EC50 = 5 * exp(ETA_EC50)
    ) |>
    transmute(
      ID = id,
      WT = wt,
      AGE = age,
      SEX = sex,
      CL, VC, Q2, V2, Q3, V3, 
      KIN, KOUT, EC50
    )
}

# simulate observations
add_obs_noise <- function(df, prop_sd_cp = 0.126, add_sd_cp = 2.09, prop_sd_cyt = 0.2, add_sd_cyt = 0.5) {
  df |> mutate(
    DV_CP  = CP * (1 + rnorm(n(), 0, prop_sd_cp)) + rnorm(n(), 0, add_sd_cp),
    DV_CYT = CYT * (1 + rnorm(n(), 0, prop_sd_cyt)) + rnorm(n(), 0, add_sd_cyt)
  )
}
