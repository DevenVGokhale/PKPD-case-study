library(tidyverse)
library(mrgsolve)

# simulate a single person ---
# compile the model
mod <- mread("./model/tgi", "true_model")

# Define dosing: 200 mg IV bolus q24h x7 doses
dose <- ev(amt = 200, ii = 24, addl = 6, cmt = 1)  # cmt = CENT

# Simulate 14 days
out <- mod %>%
  ev(dose) %>%
  mrgsim(end = 336, delta = 1) %>%
  as_tibble()

# Plot Tumor Size
out %>%
  ggplot(aes(time, TUM)) +
  geom_line(color = "firebrick", size = 1) +
  labs(title = "Tumor Size Over Time", y = "Tumor Size (mm)", x = "Time (h)") +
  theme_minimal()

# Plot Concentration
out %>%
  ggplot(aes(time, CP)) +
  geom_line(color = "steelblue", size = 1) +
  labs(title = "Plasma Concentration Over Time", y = "Concentration (mg/L)", x = "Time (h)") +
  theme_minimal()
