# requires ---
# tidyverse
# mrgsolve

# Utility functions for PKPD simulation

# Function to simulate individual characteristics
simulate_individuals <- function(n_per_arm = 10, seed = 123) {
  set.seed(seed)
  
  # Create 7 dose regimens with n_per_arm subjects each
  regimens <- tibble(regimen_id = rep(1:7, each = n_per_arm))
  
  regimens |>
    mutate(
      id = row_number(),
      age = rnorm(n(), mean = 55, sd = 12),
      age = pmax(age, 18),  # Minimum age 18
      wt = rnorm(n(), mean = 70, sd = 15),
      wt = pmax(wt, 40),    # Minimum weight 40 kg
      sex = sample(c("M", "F"), n(), replace = TRUE, prob = c(0.6, 0.4)),
      bsa = sqrt((wt * age) / 3600)  # Body surface area approximation
    )
}

# Fixed parameters
fixed_parameters <- function() {
  list(
    # PK parameters
    TVCL = 5.0,     # Typical clearance (L/h)
    TVVC = 50.0,    # Typical central volume (L)
    TVQ2 = 10.0,    # Typical inter-compartmental clearance (L/h)
    TVV2 = 100.0,   # Typical peripheral volume (L)
    TVQ3 = 2.0,     # Typical inter-compartmental clearance (L/h)
    TVV3 = 200.0,   # Typical peripheral volume (L)
    
    # PD parameters
    TVKIN = 100.0,  # Typical cytokine production rate
    TVKOUT = 0.1,   # Typical cytokine elimination rate
    TVEC50 = 10.0,  # Typical EC50 for drug effect
    
    # TGI parameters
    TVKG = 0.01,    # Typical tumor growth rate
    TVKS = 0.05,    # Typical drug killing rate
    TVLAMBDA = 0.001, # Typical resistance development rate
    TVIC50 = 5.0    # Typical IC50 for tumor killing
  )
}

# Create event schedule
create_event_schedule <- function(individuals, sim_days = 180) {
  # Define dose regimens
  regimens <- tibble(
    regimen_id = 1:7,
    dose_mg = c(25, 50, 100, 200, 400, 600, 800),
    schedule = "Q3W"  # Every 3 weeks
  )
  
  # Create dosing events for each individual
  events <- purrr::map_dfr(unique(individuals$regimen_id), function(reg_id) {
    reg_individuals <- individuals |> filter(regimen_id == reg_id)
    reg_dose <- regimens$dose_mg[regimens$regimen_id == reg_id]
    
    # Create events for each individual in this regimen
    purrr::map_dfr(reg_individuals$id, function(id) {
      # Dosing times (every 21 days)
      dose_times <- seq(0, sim_days, by = 21)
      
      # Create dosing events
      dose_events <- tibble(
        ID = id,
        time = dose_times,
        evid = 1,
        cmt = 1,
        amt = reg_dose,  # Flat dosing in mg
        rate = 0,
        addl = 0,
        ii = 0
      )
      
      # Create observation events (sparse sampling)
      obs_times <- c(0, 1, 2, 4, 8, 24, 48, 72,  # First week detailed
                     168, 336, 504, 672)  # Weekly thereafter for 4 weeks
      
      obs_events <- tibble(
        ID = id,
        time = obs_times,
        evid = 0,
        cmt = 0,
        amt = 0,
        rate = 0,
        addl = 0,
        ii = 0
      )
      
      bind_rows(dose_events, obs_events)
    })
  })
  
  events |> arrange(ID, time, desc(evid))
}

# Function to generate individual parameter data
generate_idata <- function(individuals, seed = 202) {
  set.seed(seed)
  
  # Fixed parameters
  params <- fixed_parameters()
  
  # IIV (inter-individual variability) as CV%
  iiv <- list(
    CL = 0.3,     # 30% CV
    VC = 0.25,    # 25% CV
    EC50 = 0.4,   # 40% CV
    KG = 0.25,    # 25% CV
  )
  
  n_subj <- nrow(individuals)
  
  # Generate individual parameters with log-normal distribution
  individuals |>
    mutate(
      # PK parameters
      CL = params$TVCL * exp(rnorm(n_subj, 0, iiv$CL)),
      VC = params$TVVC * exp(rnorm(n_subj, 0, iiv$VC)),
      Q2 = params$TVQ2,
      V2 = params$TVV2,
      Q3 = params$TVQ3,
      V3 = params$TVV3,
      # PD parameters
      KIN = params$TVKIN * exp(rnorm(n_subj, 0, iiv$KIN)),
      KOUT = params$TVKOUT * exp(rnorm(n_subj, 0, iiv$KOUT)),
      EC50 = params$TVEC50 * exp(rnorm(n_subj, 0, iiv$EC50)),
    ) |>
    select(ID = id, everything())  # mrgsolve expects 'ID' column
}

# Function to add observation noise using lognormal error model (proportional only)
# This matches the Stan model's lognormal likelihood
add_obs_noise <- function(df, prop_sd_cp = 0.15, prop_sd_cyt = 0.20) {
  df |>
    mutate(
      # Use lognormal error (proportional only) to match Stan model
      DV_CP = CP * exp(rnorm(n(), 0, prop_sd_cp)),
      DV_CYT = CYT * exp(rnorm(n(), 0, prop_sd_cyt)),
      # Ensure non-negative concentrations (though lognormal guarantees this)
      DV_CP = pmax(DV_CP, 0),
      DV_CYT = pmax(DV_CYT, 0)
    )
}