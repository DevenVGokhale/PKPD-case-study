## R script to simulate data using the Stan model 'three_cmt_transit_sim.stan'
## This script prepares the necessary inputs and executes the simulation.

# ------------------------------------------------------------------------------
# 1. SETUP
# ------------------------------------------------------------------------------
rm(list = ls())

# Load required libraries
library(cmdstanr)
set_cmdstan_path("~/Documents/GitHub/Torsten/cmdstan/")
library(tidyverse)

# Check if the Stan file exists
stan_file <- "./scripts_stan/three_cmt_transit_sim.stan"
if (!file.exists(stan_file)) {
  stop("Stan file not found: ", stan_file,
       ". Please ensure it's in the working directory.")
}

# Compile the Stan model
# This will create compiled C++ code for the model for faster execution later.
message("Compiling Stan model...")
mod <- cmdstan_model(stan_file)
message("Model compilation complete.")

# ------------------------------------------------------------------------------
# 2. DEFINE SIMULATION & TRIAL PARAMETERS
# (Based on pkpd_3mt_transit.R)
# ------------------------------------------------------------------------------
n_subjects <- 50
set.seed(101) # for reproducible subject weights

# Simulate subject weights
mean_wt <- 70
cv_wt <- 0.3
omega_wt <- sqrt(log(1 + cv_wt^2))
log_mean_wt <- log(mean_wt) - 0.5 * omega_wt^2
weights <- rlnorm(n_subjects, meanlog = log_mean_wt, sdlog = omega_wt)

# Define event schedule (dosing and sampling)
dose_amt <- 600      # mg, oral
dose_interval <- 72  # hours
obs_times <- seq(0, 72 * 2, by = 6) # Sampling from 0 to 144h every 6h
n_doses <- floor(max(obs_times) / dose_interval)

# Create a data frame for dosing events
doses <- tibble(
  id = 1:n_subjects,
  time = 0,
  amt = dose_amt,
  ii = dose_interval,
  addl = n_doses - 1,
  cmt = 1, # Dose goes into the first compartment (Depot)
  evid = 1,
  rate = 0,
  ss = 0
)

# Create a data frame for observation events
obs <- tibble(
  id = rep(1:n_subjects, each = length(obs_times)),
  time = rep(obs_times, times = n_subjects),
  amt = 0,
  ii = 0,
  addl = 0,
  cmt = 1, # 'cmt' for obs is a placeholder, not used by model
  evid = 0,
  rate = 0,
  ss = 0
)

# Combine and sort events to create the final schedule
events <- bind_rows(doses, obs) %>%
  arrange(id, time, desc(evid)) # Ensure dosing event comes first at t=0

# ------------------------------------------------------------------------------
# 3. PREPARE DATA LIST FOR STAN
# ------------------------------------------------------------------------------
# List of "true" parameters to simulate from
stan_data <- list(
  # Population typical values (thetas)
  CLHat = 15,
  VcHat = 60,
  Vp1Hat = 30,
  Vp2Hat = 40,
  Q2Hat = 6.0,
  Q3Hat = 5.0,
  MTTHat = 10,
  NTR = 3,
  weight_exponent_cl = 0.75,

  # IIV and Residual Error (SD values)
  omega = sqrt(log(1 + c(CL_cv = 0.3, Vc_cv = 0.2)^2)), # [omega_cl, omega_vc]
  sigma = 0.07, # Proportional error SD

  # Subject info
  nSubjects = n_subjects,
  weight = weights,
  start = which(!duplicated(events$id)),
  end = c(which(!duplicated(events$id))[-1] - 1, nrow(events)),

  # Event schedule
  nt = nrow(events),
  time = events$time,
  amt = events$amt,
  cmt = events$cmt,
  evid = events$evid,
  ii = events$ii,
  addl = events$addl,
  rate = events$rate,
  ss = events$ss,

  # Observation records
  iObsPK = which(events$evid == 0),
  nObsPK = length(which(events$evid == 0))
)


# ------------------------------------------------------------------------------
# 4. RUN SIMULATION IN STAN
# ------------------------------------------------------------------------------
message("Running simulation in Stan...")

# Use fixed_param=TRUE for simulation-only models.
# We only need 1 draw to get one simulated dataset.
sim_fit <- mod$sample(
  data = stan_data,
  fixed_param = TRUE,
  iter_warmup = 0,
  iter_sampling = 1,
  chains = 1,
  seed = 42 # for reproducible simulation
)

message("Simulation complete.")

# ------------------------------------------------------------------------------
# 5. EXTRACT AND VISUALIZE RESULTS
# ------------------------------------------------------------------------------
# Extract the simulated observations ('cObsSim')
cObsSim_vector <- as.numeric(sim_fit$draws("cObsSim"))

# Create a tidy data frame with the results
sim_data <- events %>%
  filter(evid == 0) %>% # Keep only observation records
  select(id, time) %>%
  mutate(DV = cObsSim_vector) # Add the simulated concentrations

# Plot the PK profiles for a random subset of subjects
plot_ids <- sample(1:n_subjects, 6)

pk_plot <- ggplot(sim_data %>% filter(id %in% plot_ids), aes(x = time, y = DV)) +
  geom_line(color = "dodgerblue3") +
  geom_point(color = "dodgerblue3") +
  facet_wrap(~ paste("Subject", id)) +
  labs(
    title = "Simulated PK Profiles (6 Random Subjects)",
    x = "Time (hours)",
    y = "Simulated Drug Concentration (mg/L)"
  ) +
  theme_bw()

pk_plot

# Save the plot
ggsave("simulated_pk_profiles.png", pk_plot, width = 8, height = 6)

# Optionally, save the simulated data
# write_csv(sim_data, "simulated_pk_data.csv")