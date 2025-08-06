library(cmdstanr)

mod <- cmdstan_model(stan_file = file.path(cmdstan_path(), "examples/bernoulli/bernoulli.stan"))
fit <- mod$sample(data = list(N = 10, y = rep(0:1, 5)))
print(fit)

# setting path to command torsten versions 
set_cmdstan_path("../Torsten/cmdstan/")
cmdstan_path()
cmdstan_version()

