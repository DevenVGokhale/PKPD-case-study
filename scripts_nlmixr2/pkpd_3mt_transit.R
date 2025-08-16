rm(list = ls())
library(rxode2)
library(ggplot2)
library(tidyverse)
library(MASS)
library(patchwork)
library(readr)
library(pracma)

# ------------------------------------------------------------------------------
# MODEL: 3CMT PK + Transit Absorption + Auto-Induction of CL
# ------------------------------------------------------------------------------
model <- rxode2({
  # Transit absorption
  d/dt(A0)  = -ktr * A0;
  d/dt(TR1) =  ktr * A0 - ktr * TR1;
  d/dt(TR2) =  ktr * TR1 - ktr * TR2;
  d/dt(TR3) =  ktr * TR2 - ktr * TR3;

  # 3CMT PK
  d/dt(A1) = ktr * TR3 - (CL / Vc) * A1
             - Q2 / Vc * A1 + Q2 / Vp1 * A2
             - Q3 / Vc * A1 + Q3 / Vp2 * A3;
  d/dt(A2) = Q2 / Vc * A1 - Q2 / Vp1 * A2;
  d/dt(A3) = Q3 / Vc * A1 - Q3 / Vp2 * A3;

  Cp = A1 / Vc;

})

# ------------------------------------------------------------------------------
# POPULATION PARAMETERS
# ------------------------------------------------------------------------------
theta <- c(
  CL = 15,   # L/h (baseline CL)
  Vc   = 60,      # L
  Vp1  = 30,      # L
  Vp2  = 40,      # L
  Q2   = 6.0,     # L/h
  Q3   = 5.0,     # L/h
  MTT  = 10,     # h
  NTR  = 3        # transit compartments
)

theta_wt <- 0.75
ktr <- (theta["NTR"] + 1) / theta["MTT"]

# ------------------------------------------------------------------------------
# SIMULATE SUBJECTS
# ------------------------------------------------------------------------------
n_subjects <- 50
set.seed(101)

# Simulate weights
mean_wt <- 70
cv_wt <- 0.3
omega_wt <- sqrt(log(1 + cv_wt^2))
log_mean <- log(mean_wt) - 0.5 * omega_wt^2
weights <- rlnorm(n_subjects, meanlog = log_mean, sdlog = omega_wt)
weight_df <- tibble(id = 1:n_subjects, WT = weights)

# Simulate interindividual variability
cv_iiv <- c(CL = 0.3, Vc = 0.2)
omega <- diag(log(1 + cv_iiv^2))
etas <- MASS::mvrnorm(n_subjects, mu = rep(0, length(cv_iiv)), Sigma = omega)

subject_params <- (
  weight_df |> 
    mutate(
      eta.cl = etas[, 1],
      eta.vc = etas[, 2],
      CL = theta["CL"] * (WT / 70)^theta_wt * exp(eta.cl),
      Vc   = theta["Vc"] * exp(eta.vc),
      Vp1  = theta["Vp1"],
      Vp2  = theta["Vp2"],
      Q2   = theta["Q2"],
      Q3   = theta["Q3"],
      MTT  = theta["MTT"],
      NTR  = theta["NTR"],
      ktr  = ktr,
    )
)

# ------------------------------------------------------------------------------
# EVENTS: Dosing + Sampling
# ------------------------------------------------------------------------------
dose_amt <- 600     # mg oral
dose_interval <- 72 # q24h
obs_times <- seq(0, 72*2, by = 6)  # up to 10 days
n_doses <- floor(max(obs_times)/dose_interval)

doses <- tibble(
  id = 1:n_subjects,
  time = 0,
  amt = dose_amt,
  ii = dose_interval,
  addl = n_doses - 1,
  cmt = 1,
  evid = 1
)

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
  params = subject_params |> 
            dplyr::select(id, CL, Vc, Vp1, Vp2, Q2, Q3, ktr),
  events = events,
  addDosing = TRUE
) |> 
  as_tibble() |> 
  left_join(subject_params, by = "id")

# ------------------------------------------------------------------------------
# ADD OBSERVATION NOISE
# ------------------------------------------------------------------------------
add_noise <- function(df, sd_cp = 0.07) {
  df |> mutate(
    DV_CP = Cp * (1 + rnorm(n(), 0, sd_cp)),
  )
}

sim_obs <- sim |> 
  add_noise() |> 
  dplyr::select(id, time, Cp, DV_CP, everything())

# ------------------------------------------------------------------------------
# PK METRICS
# ------------------------------------------------------------------------------
pk_metrics <- sim_obs |>
  group_by(id) |>
  summarise(
    AUC = trapz(time, DV_CP),
    Cmax = max(DV_CP),
    Tmax = time[which.max(DV_CP)[1]]
  )

# ------------------------------------------------------------------------------
# PLOT PK & PD PROFILES
# ------------------------------------------------------------------------------
if(n_subjects < 5) stop("Need more than 5 subjects!")

plot_these_ids <- sample(c(1:n_subjects), 6)

pk_plot <- ggplot(sim_obs |> filter(id %in% plot_these_ids), aes(x = time)) +
  geom_line(aes(y = DV_CP), color = "blue", alpha = 0.8) +
  geom_point(aes(y = DV_CP), color = "blue", size = 0.5) +
  facet_wrap(~id, ncol = 3, scales = "free_y") +
  labs(title = "PK Profile by Subject",
       y = "Cp (mg/L)", x = "Time (h)") +
  scale_y_continuous(transform = scales::pseudo_log_trans())+
  theme_bw()

pk_plot


pk_plot <- (
  ggplot(sim_obs |> filter(time >0), aes(x = time, group = id)) +
  geom_line(aes(y = DV_CP), color = "blue", alpha = 0.8) +  
  geom_point(aes(y = DV_CP), color = "blue", alpha = 0.5) +
  labs(title = "PK Profile (Cp) by Subject",
       y = "Cp (mg/L)", x = "Time (h)") +
  scale_y_continuous(transform = scales::pseudo_log_trans())+
  theme_bw()
)

pk_plot 

# ------------------------------------------------------------------------------
# SAVE TO FILE
# ------------------------------------------------------------------------------
write_csv(sim_obs, "data/pkpd_3cmt_transit.csv")
write_rds(sim_obs, "data/pkpd_3cmt_transit.rds")
cat("✅ Simulation complete. Saved to 'data/pkpd_3cmt_transit.{csv,rds}'\n")

