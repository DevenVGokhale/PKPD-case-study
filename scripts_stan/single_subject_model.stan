functions {
  vector pkpd_ode(real t, vector y, array[] real theta, array[] real x_r, array[] int x_i) {
    real CL   = theta[1];
    real VC   = theta[2];
    real Q2   = theta[3];
    real V2   = theta[4];
    real Q3   = theta[5];
    real V3   = theta[6];
    real KIN  = theta[7];
    real KOUT = theta[8];
    real EC50 = theta[9];

    real CP = y[1] / VC;
    real stim = CP / (EC50 + CP);

    vector[4] dydt;
    
    dydt[1] = - (CL / VC) * y[1] - (Q2 / VC) * y[1] + (Q2 / V2) * y[2]
              - (Q3 / VC) * y[1] + (Q3 / V3) * y[3];  // CENT
    dydt[2] = (Q2 / VC) * y[1] - (Q2 / V2) * y[2];     // PERIPH1
    dydt[3] = (Q3 / VC) * y[1] - (Q3 / V3) * y[3];     // PERIPH2
    dydt[4] = KIN * (1 + stim) - KOUT * y[4];          // CYT

    return dydt;
  }
}

data {
  int<lower=1> nt;                    // total events (27)
  int<lower=1> nObs;                  // number of observations (26)
  array[nObs] int<lower=1> iObs;      // observation indices
  array[nt] real time;
  array[nt] real amt;
  array[nt] int evid;
  array[nt] int cmt;
  array[nt] real rate;
  array[nt] int addl;
  array[nt] real ii;
  array[nt] int ss;
  row_vector<lower=0>[nObs] cObs;     // Combined observed concentrations (CP + CYT)
  array[nObs] int<lower=1> cmt_obs;   // Which compartment each observation is from (2=CP, 4=CYT)
  
  real rel_tol;
  real abs_tol;
  int<lower=1> max_num_steps;
}

transformed data {
  int nCmt = 4;     // 3 PK + 1 PD compartments
  int nTheta = 9;   // number of parameters
}


parameters {
  real<lower=0> CL;
  real<lower=0> VC;
  real<lower=0> Q2;
  real<lower=0> V2;
  real<lower=0> Q3;
  real<lower=0> V3;
  real<lower=0> KIN;
  real<lower=0> KOUT;
  real<lower=0> EC50;
  real<lower=0> sigma;
}

transformed parameters {
  array[nTheta] real theta = {CL, VC, Q2, V2, Q3, V3, KIN, KOUT, EC50};
  matrix[nCmt, nt] x;
  vector[nObs] cHatObs;  // Only need predictions for observations

  // Solve ODE using Torsten
  x = pmx_solve_rk45(
    pkpd_ode,
    nCmt,
    time, amt, rate, ii, evid, cmt, addl, ss,
    theta,
    rel_tol, abs_tol, max_num_steps
  );

  // Extract predictions DIRECTLY for observations using iObs
  for(i in 1:nObs) {
    if (cmt_obs[i] == 2) {          // CP observation
      cHatObs[i] = x[1, iObs[i]] / VC;  // CP = CENT/VC 
    } else if (cmt_obs[i] == 4) {   // CYT observation  
      cHatObs[i] = x[4, iObs[i]];       // CYT from compartment 4
    }
  }
}

model {
  // Informative priors from simulation setup
  CL ~ lognormal(log(0.153), 0.25);
  VC ~ lognormal(log(3.05), 0.25);
  Q2 ~ lognormal(log(0.74), 0.3);
  V2 ~ lognormal(log(1.27), 0.3);
  Q3 ~ lognormal(log(0.092), 0.3);
  V3 ~ lognormal(log(2.10), 0.3);
  KIN ~ lognormal(log(10), 0.3);
  KOUT ~ lognormal(log(1), 0.3);
  EC50 ~ lognormal(log(5), 0.3);
  sigma ~ cauchy(0, 1);

  for (i in 1:nObs) {
    if (cmt_obs[i] == 2 || cmt_obs[i] == 4) {
      if (cHatObs[i] > 0)  // Avoid log(0)
        cObs[i] ~ lognormal(log(cHatObs[i]), sigma);
    }
  }
}

generated quantities {
  array[nObs] real cObs_pred;

  for (i in 1:nObs) {
    if (cmt_obs[i] == 2 || cmt_obs[i] == 4) {
      if (cHatObs[i] > 0)
        cObs_pred[i] = lognormal_rng(log(cHatObs[i]), sigma);
      else
        cObs_pred[i] = 0;
    } else {
      cObs_pred[i] = 0;
    }
  }
}
