rm(list = ls())
library(rxode2)
library(ggplot2)
library(tidyverse)
library(MASS)
library(patchwork)

# ------------------------------------------------------------------------------
# MODEL: 3CMT IV PK + Indirect Response PD
# ------------------------------------------------------------------------------
model <- rxode2({
  d/dt(A1) = -CL/Vc * A1 
             - Q2/Vc * A1 + Q2/Vp1 * A2 
             - Q3/Vc * A1 + Q3/Vp2 * A3;
  d/dt(A2) = Q2/Vc * A1 - Q2/Vp1 * A2;
  d/dt(A3) = Q3/Vc * A1 - Q3/Vp2 * A3;

  Cp = A1 / Vc;
  
  d/dt(R) = Kin * (1 - Cp / (IC50 + Cp)) - Kout * R;
})

# ------------------------------------------------------------------------------
# POPULATION PARAMETERS (from Budha et al., 2023)
# ------------------------------------------------------------------------------
theta <- c(
  CL   = 0.0095625,   # L/h (0.153 L/day)
  Vc   = 3.05,        # L
  Vp1  = 1.27,        # L (peripheral 1)
  Vp2  = 2.10,        # L (peripheral 2)
  Q2   = 0.0308,      # L/h (0.74 L/day)
  Q3   = 0.00383,     # L/h (0.092 L/day)
  Kin  = 10,          # PD parameters (arbitrary)
  Kout = 1,
  IC50 = 2.0
)

theta_wt <- 0.75

# ------------------------------------------------------------------------------
# SIMULATE SUBJECTS
# ------------------------------------------------------------------------------
n_subjects <- 50
set.seed(101)

# Log-normal weight distribution
mean_wt <- 65  # Median weight per paper
cv_wt <- 0.3
omega_wt <- sqrt(log(1 + cv_wt^2))
log_mean <- log(mean_wt) - 0.5 * omega_wt^2
weights <- rlnorm(n_subjects, meanlog = log_mean, sdlog = omega_wt)
weight_df <- tibble(id = 1:n_subjects, WT = weights)

# IIV on CL, Vc, Kin
cv_iiv <- c(CL = 0.263, Vc = 0.167, Kin = 0.4)
omega <- diag(log(1 + cv_iiv^2))
etas <- MASS::mvrnorm(n_subjects, mu = rep(0, length(cv_iiv)), Sigma = omega)

subject_params <- weight_df |>
  mutate(
    eta.cl  = etas[, 1],
    eta.vc  = etas[, 2],
    eta.kin = etas[, 3],
    CL   = theta["CL"] * (WT / 65)^theta_wt * exp(eta.cl),
    Vc   = theta["Vc"] * exp(eta.vc),
    Vp1  = theta["Vp1"],
    Vp2  = theta["Vp2"],
    Q2   = theta["Q2"],
    Q3   = theta["Q3"],
    Kin  = theta["Kin"] * exp(eta.kin),
    Kout = theta["Kout"],
    IC50 = theta["IC50"]
  )

# ------------------------------------------------------------------------------
# DOSING EVENTS (3 mg/kg IV Q3W for 30 weeks)
# ------------------------------------------------------------------------------
dose_interval <- 21 * 24  # Q3W in hours
end_time <- 210 * 24      # 30 weeks in hours
n_doses <- floor(end_time / dose_interval)

doses <- subject_params |>
  mutate(dose_amt = 3 * WT, rate = dose_amt / 1) |>  # 1-hour infusion
  rowwise() |>
  mutate(dose_events = list(tibble(
    id = id,
    time = 0,
    amt = dose_amt,
    rate = rate,
    ii = dose_interval,
    addl = n_doses - 1,
    cmt = 1,
    evid = 1
  ))) |>
  pull(dose_events) |>
  bind_rows()

# ------------------------------------------------------------------------------
# OBSERVATION EVENTS: Dense (N=20) + Sparse (N=30) sampling
# ------------------------------------------------------------------------------

# Define IDs for dense and sparse groups
n_dense <- 20
dense_ids <- 1:n_dense
sparse_ids <- (n_dense + 1):n_subjects

# Rich sampling schedule: Early and throughout
rich_days <- c(0, 1, 2, 3, 7, 14, 21, 28, 42, 63, 84, 105, 126, 147, 168, 189, 210)
rich_times <- rich_days * 24

# Sparse sampling schedule: Only pre-dose and troughs
sparse_days <- c(0, 42, 84, 126, 168, 210)
sparse_times <- sparse_days * 24

# Build observation data
obs_rich <- expand_grid(
  id = dense_ids,
  time = rich_times
)

obs_sparse <- expand_grid(
  id = sparse_ids,
  time = sparse_times
)

# Combine observation records
obs <- bind_rows(obs_rich, obs_sparse) |>
  mutate(
    amt = 0,
    ii = 0,
    addl = 0,
    cmt = 1,
    evid = 0
  )

# Combine with dosing
events <- bind_rows(doses, obs) |> arrange(id, time, desc(evid))

# ------------------------------------------------------------------------------
# RUN SIMULATION
# ------------------------------------------------------------------------------
sim <- rxSolve(
  model,
  params = subject_params |> dplyr::select(id, CL, Vc, Vp1, Vp2, Q2, Q3, IC50, Kin, Kout),
  events = events,
  inits = c(R = 10),
  addDosing = TRUE
) |>
  as_tibble() |>
  left_join(subject_params, by = "id")

# ------------------------------------------------------------------------------
# ADD OBSERVATION NOISE
# ------------------------------------------------------------------------------
add_noise <- function(df, sd_cp = 0.07, sd_r = 0.05) {
  df |>
    mutate(
      DV_CP = Cp * (1 + rnorm(n(), 0, sd_cp)),
      DV_R  = R * (1 + rnorm(n(), 0, sd_r))
    )
}

sim_obs <- sim |>
  add_noise() |>
  dplyr::select(id, time, Cp, DV_CP, DV_R, everything())

# ------------------------------------------------------------------------------
# PLOT
# ------------------------------------------------------------------------------
plot_these_ids <- sample(c(1:n_subjects), 5)

pk_plot <- ggplot(sim_obs |> filter(id %in% plot_these_ids & time > 0), aes(x = time)) +
  geom_line(aes(y = DV_CP), color = "blue", alpha = 0.8) +
  geom_point(aes(y = DV_CP), color = "blue", size = 0.5) +
  facet_wrap(~id, ncol = 5, scales = "free_y") +
  labs(title = "PK Profile (Cp) by Subject",
       y = "Cp (mg/L)", x = "Time (h)") +
  scale_y_continuous(trans = scales::pseudo_log_trans())+
  theme_bw()

pd_plot <- ggplot(sim_obs |> filter(id %in% plot_these_ids & time >0), aes(x = time)) +
  geom_line(aes(y = DV_R), color = "darkgreen", alpha = 0.8) +
  geom_point(aes(y = DV_R), color = "darkgreen", size = 0.5) +
  facet_wrap(~id, ncol = 5, scales = "free_y") +
  labs(title = "PD Profile (Response R) by Subject",
       y = "Response (R)", x = "Time (h)") +
  theme_bw()

pk_plot / pd_plot


pk_plot <- (
  ggplot(sim_obs |> filter(time >0), aes(x = time)) +
  geom_point(aes(y = DV_CP), color = "blue", alpha = 0.5) +
  labs(title = "PK Profile (Cp) by Subject",
       y = "Cp (mg/L)", x = "Time (h)") +
  scale_y_continuous(trans = scales::pseudo_log_trans())+
  theme_bw()
)

pd_plot <- (
  ggplot(sim_obs |> filter(time >0), aes(x = time)) +
  geom_point(aes(y = DV_R), color = "darkgreen", size = 0.5) +
  labs(title = "PD Profile (Response R) by Subject",
       y = "Response (R)", x = "Time (h)") +
  theme_bw()
)

pk_plot / pd_plot

# ------------------------------------------------------------------------------
# SAVE OUTPUT
# ------------------------------------------------------------------------------
write_csv(sim_obs, "data/sbc_sim_tislelizumab.csv")
write_rds(sim_obs, "data/sbc_sim_tislelizumab.rds")
cat("✅ SBC-ready simulation complete. Saved to 'data/sbc_sim_tislelizumab.{csv,rds}'\n")
