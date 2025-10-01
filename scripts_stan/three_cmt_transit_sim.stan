// Stan model for simulating data from a three-compartment PK model
// with transit absorption. This script is intended for simulation-based
// calibration (SBC) and does not fit data. The true parameters are
// provided in the data block, and observations are generated in the
// generated quantities block.

functions {
  /**
  * ODE system for a three-compartment PK model with transit absorption.
  *
  * States (x):
  * x[1]: Amount in Depot (for transit)
  * x[2]: Amount in Transit 1
  * x[3]: Amount in Transit 2
  * x[4]: Amount in Transit 3
  * x[5]: Amount in Central Compartment
  * x[6]: Amount in Peripheral Compartment 1
  * x[7]: Amount in Peripheral Compartment 2
  *
  * Parameters (parms):
  * parms[1]: CL  (Clearance)
  * parms[2]: Vc  (Volume of Central)
  * parms[3]: Vp1 (Volume of Peripheral 1)
  * parms[4]: Vp2 (Volume of Peripheral 2)
  * parms[5]: Q2  (Inter-compartmental Clearance 1)
  * parms[6]: Q3  (Inter-compartmental Clearance 2)
  * parms[7]: ktr (Transit Rate Constant)
  */
  vector threeCptTransitModelODE(real t, vector x, array[] real parms,
                                 array[] real rdummy, array[] int idummy) {
    // Unpack parameters
    real CL  = parms[1];
    real Vc  = parms[2];
    real Vp1 = parms[3];
    real Vp2 = parms[4];
    real Q2  = parms[5];
    real Q3  = parms[6];
    real ktr = parms[7];

    // Micro-rate constants
    real k10 = CL / Vc;
    real k12 = Q2 / Vc;
    real k21 = Q2 / Vp1;
    real k13 = Q3 / Vc;
    real k31 = Q3 / Vp2;

    vector[7] dxdt;

    // Transit Absorption Compartments
    dxdt[1] = -ktr * x[1];
    dxdt[2] = ktr * x[1] - ktr * x[2];
    dxdt[3] = ktr * x[2] - ktr * x[3];
    dxdt[4] = ktr * x[3] - ktr * x[4];

    // Three-Compartment PK Model
    dxdt[5] = ktr * x[4] - (k10 + k12 + k13) * x[5] + k21 * x[6] + k31 * x[7];
    dxdt[6] = k12 * x[5] - k21 * x[6];
    dxdt[7] = k13 * x[5] - k31 * x[7];

    return dxdt;
  }
}

data {
  // Event schedule and trial design
  int<lower=1> nt;
  array[nt] real<lower=0> time;
  array[nt] real<lower=0> amt;
  array[nt] int<lower=1> cmt;
  array[nt] int<lower=0> evid;
  array[nt] real<lower=0> ii;
  array[nt] int<lower=0> addl;
  array[nt] int<lower=0> ss;
  array[nt] real rate;

  // Observation records
  int<lower=1> nObsPK;
  array[nObsPK] int<lower=1> iObsPK;

  // Subject-level information
  int<lower=1> nSubjects;
  array[nSubjects] int<lower=1> start;
  array[nSubjects] int<lower=1> end;
  array[nSubjects] real<lower=0> weight;

  // TRUE Population parameters for simulation
  real<lower=0> CLHat;
  real<lower=0> VcHat;
  real<lower=0> Vp1Hat;
  real<lower=0> Vp2Hat;
  real<lower=0> Q2Hat;
  real<lower=0> Q3Hat;
  real<lower=0> MTTHat;
  int<lower=1> NTR;
  real weight_exponent_cl;

  // TRUE IIV and residual error parameters
  vector<lower=0>[2] omega; // [omega_cl, omega_vc]
  real<lower=0> sigma;     // Proportional error SD
}

transformed data {
  int nStates = 7; // Number of states in the ODE system
  int nTheta = 7;  // Number of parameters for the ODE system
  int nIIV = 2;    // Number of parameters with IIV (CL, Vc)
  array[nSubjects] int len;

  // Length of each subject's event schedule
  for (i in 1:nSubjects) {
    len[i] = end[i] - start[i] + 1;
  }
}

parameters {
  // This block is empty because we are not fitting a model.
}

model {
  // This block is empty because we are not calculating a posterior.
}

generated quantities {
  // --- Declarations ---
  // Final simulated observations
  vector[nObsPK] cObsSim;

  // Simulated individual parameters
  array[nSubjects, nTheta] real<lower=0> parmsSim;
  matrix[nIIV, nSubjects] etaStdSim;

  // Intermediate variables for simulation
  matrix[nStates, nt] x;
  row_vector[nt] cHat;
  vector[nObsPK] cHatObs;

  // --- Simulation ---
  // 1. Simulate standard normal random effects for each subject
  for (j in 1:nIIV) {
    for (i in 1:nSubjects) {
      etaStdSim[j, i] = normal_rng(0, 1);
    }
  }

  // 2. Generate individual-specific parameters
  for (i in 1:nSubjects) {
    // Apply IIV and covariate effects
    real CL_i = CLHat * pow(weight[i] / 70.0, weight_exponent_cl)
                      * exp(omega[1] * etaStdSim[1, i]);
    real Vc_i = VcHat * exp(omega[2] * etaStdSim[2, i]);
    real ktr_i = (NTR + 1.0) / MTTHat;

    // Assemble parameter array for the ODE solver
    parmsSim[i, 1] = CL_i;
    parmsSim[i, 2] = Vc_i;
    parmsSim[i, 3] = Vp1Hat;
    parmsSim[i, 4] = Vp2Hat;
    parmsSim[i, 5] = Q2Hat;
    parmsSim[i, 6] = Q3Hat;
    parmsSim[i, 7] = ktr_i;
  }

  // 3. Solve the ODE system for all subjects
  x = pmx_solve_group_rk45(threeCptTransitModelODE, nStates, len,
                           time, amt, rate, ii, evid, cmt, addl, ss,
                           parmsSim, 1e-6, 1e-6, 10000);

  // 4. Calculate predicted concentrations
  for (i in 1:nSubjects) {
    // Concentration = Amount in Central (x[5]) / Volume of Central (parms[i,2])
    cHat[start[i]:end[i]] = x[5, start[i]:end[i]] / parmsSim[i, 2];
  }

  // 5. Extract predictions at observation times
  for (i in 1:nObsPK) {
    cHatObs[i] = cHat[iObsPK[i]];
  }

  // 6. Generate noisy observations
  for (i in 1:nObsPK) {
    // Using a lognormal error model, which corresponds to proportional
    // error on the normal scale.
    real pred = fmax(machine_precision(), cHatObs[i]);
    cObsSim[i] = exp(normal_rng(log(pred), sigma));
  }
}
