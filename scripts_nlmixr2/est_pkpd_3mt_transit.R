rm(list = ls())
library(dplyr)
library(readr)
library(tidyr)
library(nlmixr2)
library(nlmixr2est)
library(future)
library(babelmixr2)

plan(multisession, workers = parallel::detectCores() - 2)

# ------------------------------------------------------------------------------
# LOAD SIMULATED DATA
# ------------------------------------------------------------------------------
sim_obs <- read_csv("data/pkpd_3cmt_transit.csv")

# ------------------------------------------------------------------------------
# CREATE NLMIXR2 DATASET
# ------------------------------------------------------------------------------
# 1. Dosing events
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

# 2. Observation rows (Cp only)
obs_rows <- sim_obs |> 
  mutate(
    CMT = 5,  # A1 (central compartment)
    EVID = 0,
    AMT = 0,
    MDV = 0,
    DV = DV_CP
  ) |> 
  dplyr::select(ID = id, TIME = time, DV, CMT, EVID, AMT, MDV, WT)

# 3. Combine
nlmixr_data <- bind_rows(dosing_rows, obs_rows) |> 
  arrange(ID, TIME, desc(EVID))

# quick plot 
nlmixr_data |> 
  filter(EVID < 1) |> 
  ggplot(aes(x=TIME, y=DV, group=ID))+
  geom_line()+
  geom_point()+
  facet_wrap(~CMT, nrow=2, scales="free")+
  scale_y_continuous(transform = scales::pseudo_log_trans()) +
  theme_bw()

# ------------------------------------------------------------------------------
# DEFINE NLMIXR2 MODEL
# ------------------------------------------------------------------------------
pk_model <- function() {
  ini({
    tcl   <- log(14)       # Clearance
    tvc   <- log(60)       # Central Volume
    ktr   <- log((3+1)/10) # Transit Rate
    
    eta.cl ~ 0.09
    eta.vc ~ 0.04

    prop.err.cp <- 0.07
  })

  model({
    Q2   <- 6.0
    Q3   <- 5.0
    Vp1  <- 30
    Vp2  <- 40 

    CL <- exp(tcl + eta.cl) * (WT / 70)^0.75
    Vc <- exp(tvc + eta.vc)

    Cp <- A1 / Vc

    d/dt(A0)  = -ktr * A0;
    d/dt(TR1) =  ktr * A0 - ktr * TR1;
    d/dt(TR2) =  ktr * TR1 - ktr * TR2;
    d/dt(TR3) =  ktr * TR2 - ktr * TR3;

    d/dt(A1) = ktr * TR3 - (CL / Vc) * A1 
               - Q2 / Vc * A1 + Q2 / Vp1 * A2
               - Q3 / Vc * A1 + Q3 / Vp2 * A3;
    d/dt(A2) = Q2 / Vc * A1 - Q2 / Vp1 * A2;
    d/dt(A3) = Q3 / Vc * A1 - Q3 / Vp2 * A3;

    Cp ~ prop(prop.err.cp) | A1
  })
}

# ------------------------------------------------------------------------------
# FIT MODEL
# ------------------------------------------------------------------------------
fit <- nlmixr2(
  object = pk_model,
  data = nlmixr_data,
  est = "saem",
  control = saemControl(print = 50)
)

# ------------------------------------------------------------------------------
# SAVE RESULTS
# ------------------------------------------------------------------------------
print(fit)
fit$omegaR
plot(fit)
#saveRDS(fit, file = "results/fit_no_induction.rds")

#cat("\n✅ Estimation complete and saved.\n")
