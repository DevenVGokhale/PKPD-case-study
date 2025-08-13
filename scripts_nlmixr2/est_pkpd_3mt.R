rm(list = ls())
library(dplyr)
library(readr)
library(tidyr)
library(nlmixr2)
library(nlmixr2est)
library(future)
library(babelmixr2)


# Load simulated data (time in hours)
sim_obs <- read_csv("data/sbc_sim_tislelizumab.csv")


# 1. Dosing rows with correct WT
dosing_rows <- sim_obs |> 
  filter(amt > 0) |> 
  distinct(id, time, amt, .keep_all = TRUE) |> 
  mutate(
    CMT = 1,
    EVID = 1,
    AMT = amt,
    DV = NA_real_,
    MDV = 1
  ) |> 
  dplyr::select(ID = id, TIME = time, DV, CMT, EVID, AMT, MDV, WT)

# 2. Observation rows with correct CMT and WT
obs_rows <- sim_obs |> 
  pivot_longer(cols = c(DV_CP, DV_R), names_to = "endpoint", values_to = "DV") |> 
  mutate(
    CMT = case_when(
      endpoint == "DV_CP" ~ 1,   # Cp (central compartment)
      endpoint == "DV_R"  ~ 4,   # R (4th state in model)
      TRUE ~ NA_integer_  # Fail loudly
    ),
    EVID = 0,
    AMT = 0,
    MDV = 0
  ) |> 
  dplyr::select(ID = id, TIME = time, DV, CMT, EVID, AMT, MDV, WT)

# 3. Combine into final modeling dataset
nlmixr_data <- bind_rows(dosing_rows, obs_rows) |> 
  arrange(ID, TIME, desc(EVID))


# 4. Model function (rates in hours, fixed values match simulation)
tisle_sbi_model <- function() {
  ini({
    tcl    <- log(0.00956)  # CL (L/hr)
    tvc    <- log(3.05)     # Vc (L)
    #tvkin  <- log(10.0)     # Kin (/hr)

    eta.cl  ~ 0.085
    eta.vc  ~ 0.039
    #eta.kin ~ 0.144

    prop.err.cp <- 0.06
    prop.err.r  <- 0.04
  })

  model({
    # All fixed values aligned with simulation script

    q2   <- 0.0308     # L/hr
    v2   <- 1.27       # L
    q3   <- 0.00383    # L/hr
    v3   <- 2.10       # L
    kout <- 1.0        # 1/hr 
    ic50 <- 2.0        # ng/mL 

    theta_wt <- 0.75   # WT effect on CL

    cl  <- exp(tcl + eta.cl) * (WT / 65)^theta_wt
    vc  <- exp(tvc + eta.vc)
    # kin <- exp(tvkin + eta.kin)

    kin <- 10.0   
    
    cp <- A1 / vc
    inhib <- cp / (ic50 + cp)

    R(0) = 10

    d/dt(A1) = - cl / vc * A1 
               - q2 / vc * A1 + q2 / v2 * A2 
               - q3 / vc * A1 + q3 / v3 * A3
    d/dt(A2) = q2 / vc * A1 - q2 / v2 * A2
    d/dt(A3) = q3 / vc * A1 - q3 / v3 * A3
    d/dt(R)  = kin * (1 - inhib) - kout * R

    cp ~ prop(prop.err.cp) | A1
    R  ~ prop(prop.err.r)  | R
  })
}

# 5. Fit the model
plan(multisession, workers = parallel::detectCores() - 8)
fit <- nlmixr2(
  object = tisle_sbi_model,
  data = nlmixr_data,
  est = "saem",
  control = saemControl(print = 50)
)

# 6. Output results
print(fit)
# print(fit$omegaR)

vpcPlot(
  fit,
  data = nlmixr_data,
  n = 100,               # number of simulations
  stratify = "CMT",      # stratify by compartment (e.g., 1 = Cp, 4 = R)
  pred_corr = FALSE,     # or TRUE if you want prediction correction
  title = "VPC – Stratified by Compartment"
)
