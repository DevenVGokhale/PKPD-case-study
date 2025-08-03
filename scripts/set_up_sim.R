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
events_dosing <- create_event_schedule(individuals, sim_days = sim_days)


# ---- Done ----
message("✅ Simulation setup complete.")
