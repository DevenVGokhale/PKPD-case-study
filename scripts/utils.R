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
    purrr::map(as_tibble) |> 
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
    # PD parameters
    ETA_EC50 = rnorm(N, 0, sqrt(log(1 + 0.4^2)))
  )

  individuals |>
    bind_cols(eta) |>
    mutate(
      # PK params as before
      CL = 0.153 * (wt / 65)^0.565 * exp(ETA_CL),
      VC = 3.05 * (wt / 65)^0.397 * exp(ETA_VC),
      Q2 = 0.74,
      V2 = 1.27,
      Q3 = 0.092,
      V3 = 2.10,
      # PD params
      KIN = 10,
      KOUT = 1,
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

# Function to add observation noise using lognormal error model (proportional only)
# This matches the Stan model's lognormal likelihood
add_obs_noise <- function(df, prop_sd_cp = 0.15, prop_sd_cyt = 0.20) {
  df |>
    mutate(
      # Use lognormal error (proportional only) to match Stan model
      DV_CP = CP * exp(rnorm(n(), 0, prop_sd_cp)),
      DV_CYT = CYT * exp(rnorm(n(), 0, prop_sd_cyt)),
      # Ensure non-negative concentrations (though lognormal guarantees this)
      DV_CP = pmax(DV_CP, 0),
      DV_CYT = pmax(DV_CYT, 0)
    )
}

# prepare dataset to be passed on to stan/torsten
build_stan_data <- function(subject_ids=c(7,15), sim_obs_path = "data/sim_obs.rds", events_path = "data/events.rds") {
  # Load data
  sim_obs <- read_rds(sim_obs_path)
  events_raw <- read_rds(events_path)
  
  # Convert events to tibble if it's a list
  if (is.list(events_raw) && !is.data.frame(events_raw)) {
    events_data <- purrr::map_dfr(events_raw, as_tibble)
  } else {
    events_data <- events_raw
  }

  # only individuals in subject_ids will be used forward
  if (!is.null(subject_ids)) {
    subj_obs <- sim_obs |> filter(ID %in% subject_ids)
    dose_records <- events_data |> filter(ID %in% subject_ids, evid == 1) 
  } 
  
  # Get dose records (keep original ii/addl)
  dose_records <- dose_records |> 
    transmute(
      id = ID,
      time = time,
      amt = amt,
      evid = evid,
      cmt = 1L,        # Dosing compartment (like test script)
      rate = 0,
      addl = addl,     
      ii = ii,         
      ss = 0,
      DV = NA_real_,
      cmt_obs = NA_integer_
    )
  
  # CP observations - use cmt = 2 (like test script)
  cp_obs <- subj_obs |>
    filter(DV_CP > 0) |>
    transmute(
      id = ID,
      time = time,
      amt = 0,
      evid = 0L,
      cmt = 2L,        
      rate = 0,
      addl = 0,
      ii = 0,
      ss = 0,
      DV = DV_CP,
      cmt_obs = 2L 
    )
  
  # CYT observations - use cmt = 4 (for PD compartment)
  cyt_obs <- subj_obs |>
    filter(DV_CYT > 0) |>
    transmute(
      id = ID,
      time = time,
      amt = 0,
      evid = 0L,
      cmt = 4L,        # *** CHANGE: Use cmt=4 for CYT observations ***
      rate = 0,
      addl = 0,
      ii = 0,
      ss = 0,
      DV = DV_CYT,
      cmt_obs = 4L
    )
  
  # Rest stays the same...
  all_events <- bind_rows(dose_records, cp_obs, cyt_obs) |>
    arrange(id, time, desc(evid))

  start_end_data <- (
    all_events |> 
      mutate(row_index = row_number()) |> 
      group_by(id) |> 
      summarise(
        start = min(row_index),
        end = max(row_index)) |> 
      ungroup()  
  )

  iObs <- which(all_events$evid == 0)
  nObs <- length(iObs)
  nt <- nrow(all_events)
  
  cObs <- all_events$DV[iObs]
  cmt_obs <- all_events$cmt_obs[iObs]
  
  stan_data <- list(
    nSubjects = nrow(start_end_data),
    nt = nt,
    nObs = nObs,
    iObs = iObs,
    time = all_events$time,
    cObs = cObs,
    cmt_obs = cmt_obs,
    start = start_end_data$start,
    end = start_end_data$end,
    WT = unique(subj_obs$WT),
    amt = all_events$amt,
    rate = all_events$rate,
    cmt = all_events$cmt,      # Now will be [1,2,4,2,4,2,4,...]
    evid = all_events$evid,
    ii = as.integer(all_events$ii),
    addl = as.integer(all_events$addl),
    ss = all_events$ss,
    rel_tol = 1e-6,
    abs_tol = 1e-6,
    max_num_steps = 5000
  )
  
  return(stan_data)
}

# pop model 
build_stan_data_pop <- function(subject_ids = c(7, 15),
                                sim_obs_path = "data/sim_obs.rds",
                                events_path = "data/events.rds") {
  library(dplyr)
  library(readr)
  library(purrr)
  
  # Load data
  sim_obs <- read_rds(sim_obs_path)
  events_raw <- read_rds(events_path)
  
  # Convert events to tibble if it's a list
  if (is.list(events_raw) && !is.data.frame(events_raw)) {
    events_data <- purrr::map_dfr(events_raw, as_tibble)
  } else {
    events_data <- events_raw
  }

  # Subset only subjects of interest
  subj_obs <- sim_obs |> filter(ID %in% subject_ids)
  dose_records <- events_data |> filter(ID %in% subject_ids, evid == 1)

  # Dose rows
  dose_records <- dose_records |> 
    transmute(
      id = ID,
      time, amt, evid,
      cmt = 1L,
      rate = 0,
      addl, ii,
      ss = 0,
      DV = NA_real_,
      cmt_obs = NA_integer_
    )

  # Observation rows
  cp_obs <- subj_obs |>
    filter(DV_CP > 0) |>
    transmute(
      id = ID,
      time, amt = 0, evid = 0L,
      cmt = 2L, rate = 0, addl = 0L, ii = 0, ss = 0,
      DV = DV_CP, cmt_obs = 2L
    )

  cyt_obs <- subj_obs |>
    filter(DV_CYT > 0) |>
    transmute(
      id = ID,
      time, amt = 0, evid = 0L,
      cmt = 4L, rate = 0, addl = 0L, ii = 0, ss = 0,
      DV = DV_CYT, cmt_obs = 4L
    )

  # Combine and sort
  all_events <- bind_rows(dose_records, cp_obs, cyt_obs) |> 
    arrange(id, time, desc(evid))

  # Compute indices for solver
  start_end_data <- all_events |>
    mutate(row_index = row_number()) |>
    group_by(id) |>
    summarise(start = min(row_index), end = max(row_index), .groups = "drop")

  # Observation indices
  pk_idx <- which(all_events$cmt_obs == 2)
  pd_idx <- which(all_events$cmt_obs == 4)

  stan_data <- list(
    nSubjects = nrow(start_end_data),
    nt = nrow(all_events),
    nObsPK = length(pk_idx),
    nObsPD = length(pd_idx),
    iObsPK = pk_idx,
    iObsPD = pd_idx,
    time = all_events$time,
    amt = all_events$amt,
    rate = all_events$rate,
    ii = as.numeric(all_events$ii),  # real
    evid = all_events$evid,
    cmt = all_events$cmt,
    addl = all_events$addl,
    ss = all_events$ss,
    cObs = all_events$DV[pk_idx],
    pdObs = all_events$DV[pd_idx],
    # cmt_obs = all_events$cmt_obs,
    start = start_end_data$start,
    end = start_end_data$end,
    WT = subj_obs |> group_by(ID) |> summarise(WT = first(WT)) |> pull(WT),
    rel_tol = 1e-6,
    abs_tol = 1e-6,
    max_num_steps = 5000
  )

  return(stan_data)
}
