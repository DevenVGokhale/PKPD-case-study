rm(list = ls())
library(rxode2)
library(ggplot2)
library(tidyverse)
library(purrr)
library(MASS)
library(readr)
library(patchwork)

# ------------------------------------------------------------------------------
# MODEL DEFINITION (RxODE tutorial style: indirect response inhibition)
# ------------------------------------------------------------------------------
model <- rxode2({
  d/dt(A1) = -CL / Vc * A1;
  Cp = A1 / Vc;
  d/dt(R) = Kin * (1 - Cp / (IC50 + Cp)) - Kout * R;
})

# ------------------------------------------------------------------------------
# POPULATION PARAMETERS + IIV
# ------------------------------------------------------------------------------
theta_wt <- 0.75

theta <- c(
  CL = 1.5,
  Vc = 15,
  Kin = 10,
  Kout = 1,
  IC50 = 2.0
)

n_subjects <- 50
set.seed(101)

# ---- Simulating weights ----
# Assume log-normal with mean = 70 kg and CV = 30%
mean_wt <- 70
cv_wt <- 0.3
omega_wt <- sqrt(log(1 + cv_wt^2))
log_mean <- log(mean_wt) - 0.5 * omega_wt^2  # Adjusted log-mean for unbiased back-transformation
weights <- rlnorm(n_subjects, meanlog = log_mean, sdlog = omega_wt)
weight_df <- tibble(id = 1:n_subjects, WT = weights)

# ---- simulating iivs on the parameters ----
cv_iiv <- c(CL = 0.3, Vc = 0.2, Kout = 0.4)
omega <- diag(log(1 + cv_iiv^2))

etas <- MASS::mvrnorm(n_subjects, mu = rep(0, length(cv_iiv)), Sigma = omega)

subject_params <- (
  weight_df |> 
    mutate(
      eta.cl = etas[, 1],
      eta.vc = etas[, 2],
      eta.kin = etas[, 3],
      CL   = theta["CL"] * (WT / 70)^theta_wt * exp(eta.cl),
      Vc   = theta["Vc"] * exp(eta.vc),
      IC50 = theta["IC50"],
      Kin  = theta["Kin"] * exp(eta.kin),
      Kout = theta["Kout"]
    )
)
  

# ------------------------------------------------------------------------------
# DOSING + OBSERVATION SCHEDULE
# ------------------------------------------------------------------------------
dose_amt <- 100     # mg
dose_interval <- 48 # q48h
obs_times <- seq(0, 240, by = 4) # hourly sampling over 10 days
n_doses <- floor(max(obs_times)/dose_interval)        # total doses

# Dosing events
doses <- tibble(
  id = 1:n_subjects,
  time = 0,
  amt = dose_amt,
  ii = dose_interval,
  addl = n_doses - 1,
  cmt = 1,
  evid = 1
)

# Observation events
obs <- tibble(
  id = rep(1:n_subjects, each = length(obs_times)),
  time = rep(obs_times, times = n_subjects),
  amt = 0,
  ii = 0,
  addl = 0,
  cmt = 1,
  evid = 0
)

events <- bind_rows(doses, obs) |> arrange(id, time, desc(evid))

# ------------------------------------------------------------------------------
# RUN SIMULATION
# ------------------------------------------------------------------------------
sim <- rxSolve(
  model,
  params = subject_params |> select(id, CL, Vc, IC50, Kin, Kout),
  events = events,
  inits = c(R = 10),  # explicit initial conditions
  addDosing = TRUE
) |> 
  as_tibble() |> 
  left_join(subject_params, by = "id")

# ------------------------------------------------------------------------------
# ADD OBSERVATION NOISE (proportional error)
# ------------------------------------------------------------------------------
add_noise <- function(df, sd_cp = 0.07, sd_r = 0.05) {
  df |> mutate(
    DV_CP = Cp * (1 + rnorm(n(), 0, sd_cp)),
    DV_R = R * (1+ rnorm(n(), 0, sd_r))
  )
}

sim_obs <- sim |> 
  add_noise() |> 
  dplyr::select(id, time, Cp, DV_CP, DV_R, everything())



# ------------------------------------------------------------------------------
# PLOT: PK and PD profiles grouped by subject
# ------------------------------------------------------------------------------

if(n_subjects < 5) {
  stop("Need more than 5 subjects!")
}

plot_these_ids <- sample(c(1:n_subjects), 5)

pk_plot <- ggplot(sim_obs |> filter(id %in% plot_these_ids), aes(x = time)) +
  geom_line(aes(y = DV_CP), color = "blue", alpha = 0.8) +
  geom_point(aes(y = DV_CP), color = "blue", size = 0.5) +
  facet_wrap(~id, ncol = 5, scales = "free_y") +
  labs(title = "PK Profile (Cp) by Subject",
       y = "Cp (mg/L)", x = "Time (h)") +
  theme_bw()

pd_plot <- ggplot(sim_obs |> filter(id %in% plot_these_ids), aes(x = time)) +
  geom_line(aes(y = DV_R), color = "darkgreen", alpha = 0.8) +
  geom_point(aes(y = DV_R), color = "darkgreen", size = 0.5) +
  facet_wrap(~id, ncol = 5, scales = "free_y") +
  labs(title = "PD Profile (Response R) by Subject",
       y = "Response (R)", x = "Time (h)") +
  theme_bw()

pk_plot / pd_plot


# ------------------------------------------------------------------------------
# SAVE OUTPUT
# ------------------------------------------------------------------------------
write_csv(sim_obs, "data/pkpd_indirect_inhibition.csv")
write_rds(sim_obs, "data/pkpd_indirect_inhibition.rds")

cat("✅ Simulation complete. Saved to 'data/pkpd_indirect_inhibition.{csv,rds}'\n")
