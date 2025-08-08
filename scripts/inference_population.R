rm(list = ls())
library(tidyverse)
library(cmdstanr)
library(posterior)

source("./scripts/utils.R")
#debugonce(build_stan_data)
stan_data <- build_stan_data_pop(subject_ids=c(21, 15))

# set the path to the model 
pop_model_path <- "scripts_stan/population_model.stan"

# Compile the model
pop_mod <- cmdstan_model(pop_model_path)

fit <- pop_mod$sample(
  data = stan_data,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  refresh = 100,
  seed = 123
)

# Check the fit
print(fit$summary())
