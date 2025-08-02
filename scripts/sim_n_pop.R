rm(list = ls())
library(tidyverse)
library(mrgsolve)

source("./scripts/utils.R")
source("./scripts/set_up_sim.R")

# simulate a a pop of size n ---
# compile the model
mod <- mread("./model/true_model")

out <- mod |> mrgsim_d(data = events, idata = idata) |> as_tibble()
