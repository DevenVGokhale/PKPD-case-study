# 1-Compartment PK-PD Simulation with No IIV Correlation
# Uses rxode2 for simulation

library(rxode2)
library(dplyr)
library(ggplot2)
library(mvtnorm)

# Define the PK-PD model
pkpd_model <- rxode2({
  # PK parameters with individual variability
  CLi = CL * exp(eta_CL)
  Vi = VC * exp(eta_VC)
  
  # PD parameters with individual variability  
  EC50i = EC50 * exp(eta_EC50)
  
  # PK model (1-compartment IV)
  d/dt(CENT) = -CLi/Vi * CENT
  CP = CENT / Vi
  
  # PD model (Emax model)
  EFFECT = EMAX * CP / (EC50i + CP)
  
  # Add residual error
  CP_obs = CP * (1 + eps_CP)
  EFFECT_obs = EFFECT * (1 + eps_PD)
})

# Population parameters (realistic values)
theta <- c(
  CL = 10,      # Clearance (L/h) 
  VC = 50,      # Volume of distribution (L)
  EMAX = 100,   # Maximum effect
  EC50 = 5      # Concentration for 50% effect (mg/L)
)

# IIV parameters (CV%) - NO CORRELATION
omega_no_corr <- c(
  eta_CL = 0.3,     # 30% CV on CL
  eta_VC = 0.2,     # 20% CV on VC  
  eta_EC50 = 0.4    # 40% CV on EC50
)

# Residual error
sigma <- c(
  eps_CP = 0.15,    # 15% proportional error on CP
  eps_PD = 0.10     # 10% proportional error on PD
)

# Create individual parameters (10 subjects, no correlation)
set.seed(123)
n_subjects <- 10

# Generate uncorrelated random effects
eta_matrix_no_corr <- matrix(0, nrow = n_subjects, ncol = 3)
eta_matrix_no_corr[,1] <- rnorm(n_subjects, 0, omega_no_corr[1])  # eta_CL
eta_matrix_no_corr[,2] <- rnorm(n_subjects, 0, omega_no_corr[2])  # eta_VC
eta_matrix_no_corr[,3] <- rnorm(n_subjects, 0, omega_no_corr[3])  # eta_EC50

# Create event table (IV infusion dosing)
# 3 doses: 100 mg at 0h, 200 mg at 24h, 150 mg at 48h
# Each dose infused over 1 hour

ev <- et() %>%
  et(amt = 100, time = 0, rate = 100, cmt = "CENT") %>%    # 100mg over 1h
  et(amt = 200, time = 24, rate = 200, cmt = "CENT") %>%   # 200mg over 1h  
  et(amt = 150, time = 48, rate = 150, cmt = "CENT") %>%   # 150mg over 1h
  et(seq(0, 72, by = 1))  # Observations every hour for 72h

# Create individual parameter data frame
idata_no_corr <- data.frame(
  id = 1:n_subjects,
  CL = theta["CL"],
  VC = theta["VC"], 
  EMAX = theta["EMAX"],
  EC50 = theta["EC50"],
  eta_CL = eta_matrix_no_corr[,1],
  eta_VC = eta_matrix_no_corr[,2],
  eta_EC50 = eta_matrix_no_corr[,3],
  eps_CP = sigma["eps_CP"],
  eps_PD = sigma["eps_PD"]
)

# Run simulation
sim_no_corr <- rxSolve(pkpd_model, ev, idata_no_corr, addDosing = TRUE)

# Convert to data frame and add true individual parameters
sim_df_no_corr <- as.data.frame(sim_no_corr) %>%
  mutate(
    scenario = "No Correlation",
    true_CL = theta["CL"] * exp(eta_CL),
    true_VC = theta["VC"] * exp(eta_VC),
    true_EC50 = theta["EC50"] * exp(eta_EC50)
  )

# Save simulation results
write.csv(sim_df_no_corr, "data/pkpd_sim_no_correlation.csv", row.names = FALSE)

# Summary of individual parameters
cat("=== PK-PD Simulation (No IIV Correlation) ===\n")
cat("Population parameters:\n")
print(theta)
cat("\nIIV (CV%):\n")
print(omega_no_corr)
cat("\nResidual error:\n")
print(sigma)

# Individual parameter summary
individual_params_no_corr <- idata_no_corr %>%
  mutate(
    CL_individual = CL * exp(eta_CL),
    VC_individual = VC * exp(eta_VC),
    EC50_individual = EC50 * exp(eta_EC50)
  ) %>%
  select(id, CL_individual, VC_individual, EC50_individual)

cat("\nIndividual parameters:\n")
print(individual_params_no_corr)

# Check correlations between individual parameters
cat("\nCorrelation matrix (individual parameters):\n")
cor_matrix_no_corr <- cor(individual_params_no_corr[,-1])
print(round(cor_matrix_no_corr, 3))

cat("\nSimulation completed! Data saved to 'data/pkpd_sim_no_correlation.csv'\n")