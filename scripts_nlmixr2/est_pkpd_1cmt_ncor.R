rm(list = ls())
library(dplyr)
library(readr)
library(tidyr)
library(nlmixr2)
library(nlmixr2est)
library(future)
library(babelmixr2)

plan(multisession, workers = parallel::detectCores() - 8)

# Load simulated data
sim_obs <- read_csv("data/pkpd_indirect_inhibition.csv")

# 1. Extract dosing events
dosing_rows <- sim_obs |> 
  filter(amt > 0) |> 
  distinct(id, time, amt, .keep_all = TRUE) |> 
  mutate(
    CMT = 1,    # Dose goes into A1
    EVID = 1,
    AMT = amt,
    DV = NA_real_,
    MDV = 1,
    WT = NA_real_
  ) |> 
  dplyr::select(ID = id, TIME = time, DV, CMT, EVID, AMT, MDV, WT)

# 2. Reshape observation rows (DV_CP → Cp, DV_R → R)
obs_rows <- sim_obs |> 
  pivot_longer(cols = c(DV_CP, DV_R), names_to = "endpoint", values_to = "DV") |>
  mutate(
    CMT = case_when(
      endpoint == "DV_CP" ~ 1,  # A1/Cp observation
      endpoint == "DV_R" ~ 2,   # R observation
    ),
    EVID = 0,
    AMT = 0,
    MDV = 0
  ) |> 
  dplyr::select(ID = id, TIME = time, DV, CMT, EVID, AMT, MDV, WT)

# 3. Combine
nlmixr_data <- bind_rows(dosing_rows, obs_rows) |> 
  arrange(ID, TIME, desc(EVID))

pkpd_model <- function() {
  ini({
    tcl    <- log(1.5)   # log(CL)
    tvc    <- log(15)    # log(Vc)
    tvkin  <- log(10.0)   # log(KIN)
    
    theta_wt <- 0.75      # Effect of WT on CL

    eta.cl  ~ 0.085        # ≈ 30% CV
    eta.vc  ~ 0.039       # ≈ 20% CV
    eta.kin ~ 0.144       # ≈ 40% CV

    prop.err.cp <- 0.06  # SD for Cp
    prop.err.r  <- 0.04  # SD for R
  })

  model({
    cl   <- exp(tcl + eta.cl) * (WT / 70)^theta_wt
    vc   <- exp(tvc + eta.vc)
    ic50 <- 2

    kin  <- exp(tvkin + eta.kin)
    kout <- 1
    

    cp <- A1 / vc
    pd <- 1 - cp / (ic50 + cp)
    
    R(0) = 10
    d/dt(A1) = -cl * cp
    d/dt(R)  = kin * pd - kout * R

    cp ~ prop(prop.err.cp) | A1
    R  ~ prop(prop.err.r)  | R
  })
}


fit <- nlmixr2(
  object = pkpd_model,
  data = nlmixr_data,
  est = "saem",
  control = saemControl(print = 50)
)

fit

fit$omegaR

# plot(fit)


# library(xpose.nlmixr2)
# xpdb <- xpose_data_nlmixr(fit)
# dv_vs_pred(xpdb)
# res_vs_pred(xpdb)

# fit$eta |> 
#   as_tibble(rownames = "ID") |> 
#   pivot_longer(cols = -ID, names_to = "parameter", values_to = "eta") |> 
#   ggplot(aes(x = eta)) +
#   geom_histogram(bins = 20, fill = "steelblue", color = "white", alpha = 0.8) +
#   geom_density(color = "darkred", linewidth = 1.2) +
#   facet_wrap(~parameter, scales = "free") +
#   theme_bw() +
#   labs(
#     title = "Distribution of Empirical Bayes Estimates (EBEs)",
#     x = "ETA value",
#     y = "Count"
#   )
