install.packages(c("nlmixr2", "nlmixr2est", "babelmixr2"), dependencies = TRUE)

library(nlmixr2)       # Load THIS first
library(nlmixr2est)    # Estimation backend
library(babelmixr2)    # Backend glue

# Define a minimal model (no data needed)
model_test <- function() {
  ini({
    tCL <- log(1)
    eta.CL ~ 0.1
    prop.err <- 0.1
  })

  model({
    CL <- exp(tCL + eta.CL)
    d/dt(A1) = -CL * A1
    if (CMT == 1) {
      F <- Cp
      Y <- F * (1 + prop.err)
    }
    
  })
}

print(model_test())
