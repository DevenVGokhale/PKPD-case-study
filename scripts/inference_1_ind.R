rm(list = ls())
library(tidyverse)
library(cmdstanr)
library(posterior)

library(tidyverse)

build_stan_data <- function(subject_id, sim_obs_path = "data/sim_obs.rds", events_path = "data/events.rds") {
  # Load data
  sim_obs <- read_rds(sim_obs_path)
  events_data <- read_rds(events_path)

  # Filter subject's data
  subj_obs <- sim_obs |> filter(ID == subject_id)
  
  # Get dose records (keep original ii/addl)
  dose_records <- events_data |> 
    filter(ID == subject_id, evid == 1) |>
    transmute(
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
      time = time,
      amt = 0,
      evid = 0L,
      cmt = 2L,        # *** CHANGE: Use cmt=2 like test script ***
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
    arrange(time, desc(evid))
  
  iObs <- which(all_events$evid == 0)
  nObs <- length(iObs)
  nt <- nrow(all_events)
  
  cObs <- all_events$DV[iObs]
  cmt_obs <- all_events$cmt_obs[iObs]
  
  stan_data <- list(
    nt = nt,
    nObs = nObs,
    iObs = iObs,
    time = all_events$time,
    cObs = cObs,
    cmt_obs = cmt_obs,
    amt = all_events$amt,
    rate = all_events$rate,
    cmt = all_events$cmt,      # Now will be [1,2,4,2,4,2,4,...]
    evid = all_events$evid,
    ii = all_events$ii,
    addl = all_events$addl,
    ss = all_events$ss,
    rel_tol = 1e-6,
    abs_tol = 1e-6,
    max_num_steps = 1000
  )
  
  return(stan_data)
}

subject_id <- 7

stan_data <- build_stan_data(subject_id)

# set the path to the model 
model_path <- "scripts_stan/single_subject_model.stan"

# Compile the model
mod <- cmdstan_model(model_path)

# fitting the model to the data 
fit <- mod$sample(
  data = stan_data,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  refresh = 100,
  seed = 123
)

# Check the fit
print(fit$summary())
