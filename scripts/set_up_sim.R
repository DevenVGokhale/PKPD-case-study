#!/usr/bin/env Rscript
# ---- Depends on ----
# libraries ----
# tidyverse
# mrgsolve

# scripts ----
# utils.R 

# ---- Parameters ----
n_per_regimen <- 10
sim_days <- 180
set.seed(42)

# ---- Step 1: Simulate Individuals ----
individuals <- simulate_individuals(n_per_arm = n_per_regimen)

# ---- Step 2: Generate idata ----
idata <- generate_idata(individuals)

# ---- Step 3: Create event schedule ----
events <- create_event_schedule(individuals, sim_days = sim_days)

# ---- Step 4: Save inputs for mrgsolve ----
write_rds(individuals, file = "data/individuals.rds")
write_rds(idata, file = "data/idata.rds")
write_rds(events, file = "data/events.rds")

# ---- Done ----
message("✅ Simulation setup complete. Saved to 'data/'.")
