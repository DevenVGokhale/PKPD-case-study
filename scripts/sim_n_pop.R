rm(list = ls())

library(tidyverse)
library(mrgsolve)
source("./scripts/utils.R")  
source("./scripts/set_up_sim.R")

# ---- 1. Set up population ----
n_per_regimen <- 10
sim_days <- 180

# ---- 2. Add realistic sparse observation times ----
sampling_times <- tibble(
  time = c(0.5, 1, 2, 4, 8, 24, 48, 72, 120, 168, 336, 504, 672)  # q3w sampling
)

obs_events <- individuals |> 
  crossing(sampling_times) |> 
  mutate(
    amt = 0, evid = 0, cmt = 0, mdv = 0
  ) |> 
  select(ID = id, time, amt, evid, cmt, mdv)

# mrgsolve expects all event data to be a single data frame
events <- bind_rows(events_dosing, obs_events) |> arrange(ID, time, desc(evid))

# ---- 3. Simulate model ----
mod <- mread("true_model", project = "model")
sim_dense <- (
  mod |> 
    mrgsim_d(data = events, idata = idata, 
             carry_out = c("evid", "amt", "cmt")) |> 
    as_tibble())


# renaming individuals
individuals <- (
  individuals |> 
    rename(ID = id, REGIMEN = regimen, SEX=sex, WT=wt, AGE=age)
)

sim_obs <- sim_dense |>
  filter(evid == 0) |>  # only observations
  select(ID, time, evid, cmt, amt, CP, CYT) |>
  left_join(individuals, by = "ID") |>
  left_join(idata |> select(-WT, -AGE, -SEX), by = "ID") |>  # avoid duplication
  add_obs_noise() |>
  mutate(
    DV = DV_CP,     # Choose which DV you're modeling (CP vs CYT)
    MDV = 0,        # Missing DV flag for NONMEM
    EVID = evid,    # duplicate for clarity
    CMT = 2         # assumes CP is in compartment 2 (can adjust)
  ) |>
  select(ID, time, DV, EVID, AMT = amt, CMT, MDV, everything())

# ---- 5. Save output ----
write_rds(individuals, "data/individuals.rds")
write_rds(idata, "data/idata.rds")
write_rds(events, "data/events.rds")
write_rds(sim_dense, "data/sim_dense.rds")
write_rds(sim_obs, "data/sim_obs.rds")
write_csv(sim_obs, "data/sim_obs.csv")  # optional for inspection

message("✅ Synthetic data generation complete.")
