library(tidyverse)
library(mrgsolve)
set.seed(42)

# simulate a a pop of size n ---
# compile the model
mod <- mread("./model/tgi", "true_model")

# Define variability
eta_vals <- function(n) {
  tibble(
    ID = 1:n,
    ETA_CL = rnorm(n, 0, 0.2),
    ETA_V1 = rnorm(n, 0, 0.2),
    ETA_KG = rnorm(n, 0, 0.3),
    ETA_KS = rnorm(n, 0, 0.3)
  )
}

# Define dose regimens
dose_arm_a <- ev(amt = 100, ii = 24, addl = 6, cmt = 1)
dose_arm_b <- ev(amt = 200, ii = 24, addl = 6, cmt = 1)
dose_arm_c <- ev(amt = 300, ii = 48, addl = 3, cmt = 1)

# Assign arms
arms <- tribble(
  ~ARM, ~DOSE,
  "A", dose_arm_a,
  "B", dose_arm_b,
  "C", dose_arm_c
)

n_patients <- 30

# Combine all arms into one simulation
sim_all <- purrr::map_dfr(1:nrow(arms), function(i) {
  arm <- arms$ARM[i]
  dose <- arms$DOSE[[i]]
  
  etas <- eta_vals(n_patients) %>% mutate(ARM = arm)
  
  sim <- mod %>%
    idata_set(etas) %>%   # <- key change here
    ev(dose) %>%
    mrgsim(end = 336, delta = 12) %>%
    as_tibble()
  
  left_join(sim, etas, by = "ID")
})

sim_all %>%
  ggplot(aes(time, TUM, group = ID)) +
  geom_line(alpha = 0.2) +
  facet_wrap(~ARM) +
  labs(title = "Tumor Response by Dosing Arm", y = "Tumor Size (mm)", x = "Time (h)") +
  theme_minimal()

sim_all %>%
  group_by(ARM, time) %>%
  summarise(
    median = median(TUM),
    lwr = quantile(TUM, 0.1),
    upr = quantile(TUM, 0.9),
    .groups = "drop"
  ) %>%
  ggplot(aes(time, median, fill = ARM, color = ARM)) +
  geom_line() +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2, color = NA) +
  facet_wrap(~ARM) +
  labs(title = "Tumor Size (Median ± 10th–90th Percentiles)", y = "Tumor Size", x = "Time") +
  theme_minimal()
