rm(list = ls())
library(tidyverse)
library(mrgsolve)

source("./scripts/utils.R")
source("./scripts/set_up_sim.R")
source("./scripts/sim_n_pop.R")

theme_set(theme_bw())

plt1 <- (
  sim_all |> 
    ggplot(aes(x = time, y = CP_obs, group = ID)) +
    geom_line() +
    #facet_wrap(~SEX, scales = "free_y") +
    labs(title = "Drug Concentration (CP) by Individual", 
        y = "CP (µg/mL)", x = "Time (days)")
      )
