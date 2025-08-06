functions {
  vector pkpd_ode(real t, vector y, array[] real parms,
                 array[] real x_r, array[] int x_i) {
    real CL = parms[1];
    real VC = parms[2];
    real Q2 = parms[3];
    real V2 = parms[4];
    real Q3 = parms[5];
    real V3 = parms[6];
    real KIN = parms[7];
    real KOUT = parms[8];
    real EC50 = parms[9];

    real CP = y[1] / VC;
    real stim = CP / (EC50 + CP);

    vector[4] dydt;
    dydt[1] = - (CL / VC) * y[1] - (Q2 / VC) * y[1] + (Q2 / V2) * y[2]
              - (Q3 / VC) * y[1] + (Q3 / V3) * y[3];
    dydt[2] = (Q2 / VC) * y[1] - (Q2 / V2) * y[2];
    dydt[3] = (Q3 / VC) * y[1] - (Q3 / V3) * y[3];
    dydt[4] = KIN * (1 + stim) - KOUT * y[4];

    return dydt;
  }
}

data {
  int<lower=1> nt;
  int<lower=1> nObsPK;
  int<lower=1> nObsPD;
  array[nObsPK] int<lower=1> iObsPK;
  array[nObsPD] int<lower=1> iObsPD;
  array[nt] real amt;
  array[nt] int cmt;
  array[nt] int evid;
  array[nt] real time;
  array[nt] real ii;
  array[nt] int addl;
  array[nt] int ss;
  array[nt] real rate;
  vector<lower=0>[nObsPK] cObs;
  vector<lower=0>[nObsPD] pdObs;

  int<lower=1> nSubjects;
  array[nSubjects] int<lower=1> start;
  array[nSubjects] int<lower=1> end;
  array[nSubjects] real<lower=0> WT;

  real rel_tol;
  real abs_tol;
  int<lower=1> max_num_steps;
}

transformed data {
  int nTheta = 9;
  int nIIV = 3;
  int nCmt = 4;
  array[nSubjects] int len;
  for (i in 1:nSubjects)
    len[i] = end[i] - start[i] + 1;
}

parameters {
  real<lower=0> CLHat;
  real<lower=0> VCHat;
  real<lower=0> Q2;
  real<lower=0> V2;
  real<lower=0> Q3;
  real<lower=0> V3;
  real<lower=0> KIN;
  real<lower=0> KOUT;
  real<lower=0> EC50Hat;
  real<lower=0> sigmaPK;
  real<lower=0> sigmaPD;

  real<lower=0, upper=1> theta_CL_WT;
  real<lower=0, upper=1> theta_VC_WT;

  cholesky_factor_corr[nIIV] L;
  vector<lower=0>[nIIV] omega;
  matrix[nIIV, nSubjects] etaStd;
}

transformed parameters {
  array[nSubjects, nTheta] real<lower = 0> theta;
  matrix[nCmt, nt] x;
  vector[nObsPK] cHatObs;
  vector[nObsPD] pdHatObs;

  vector[nIIV] thetaHat = [log(CLHat), log(VCHat), log(EC50Hat)]';
  matrix[nSubjects, nIIV] thetaM =
    exp(rep_matrix(thetaHat, nSubjects)' + (diag_pre_multiply(omega, L) * etaStd)');

  for (i in 1:nSubjects) {
    theta[i, 1] = thetaM[i, 1] * pow(WT[i] / 65, theta_CL_WT);
    theta[i, 2] = thetaM[i, 2] * pow(WT[i] / 65, theta_VC_WT);
    theta[i, 3] = Q2;
    theta[i, 4] = V2;
    theta[i, 5] = Q3;
    theta[i, 6] = V3;
    theta[i, 7] = KIN;
    theta[i, 8] = KOUT;
    theta[i, 9] = thetaM[i, 3];
  }

  x = pmx_solve_group_rk45(
    pkpd_ode, nCmt,
    len,
    time, amt, rate, ii,
    evid, cmt, addl, ss,
    theta,
    rel_tol, abs_tol, max_num_steps
  );

  for (i in 1:nObsPK)
    cHatObs[i] = x[1, iObsPK[i]] / theta[1, 2];
  for (i in 1:nObsPD)
    pdHatObs[i] = x[4, iObsPD[i]];
}

model {
  CLHat ~ lognormal(log(0.153), 0.25);
  VCHat ~ lognormal(log(3.05), 0.25);
  Q2 ~ lognormal(log(0.74), 0.3);
  V2 ~ lognormal(log(1.27), 0.3);
  Q3 ~ lognormal(log(0.092), 0.3);
  V3 ~ lognormal(log(2.10), 0.3);
  KIN ~ lognormal(log(10), 0.3);
  KOUT ~ lognormal(log(1), 0.3);
  EC50Hat ~ lognormal(log(5), 0.3);
  sigmaPK ~ cauchy(0, 1);
  sigmaPD ~ cauchy(0, 1);

  theta_CL_WT ~ normal(0.75, 0.2);
  theta_VC_WT ~ normal(1.0, 0.2);

  L ~ lkj_corr_cholesky(1);
  to_vector(etaStd) ~ normal(0, 1);
  omega ~ cauchy(0, 0.5);

  cObs ~ lognormal(log(cHatObs), sigmaPK);
  pdObs ~ lognormal(log(pdHatObs), sigmaPD);
}

generated quantities {
  vector[nObsPK] cObs_pred;
  vector[nObsPD] pdObs_pred;

  for (i in 1:nObsPK)
    cObs_pred[i] = lognormal_rng(log(cHatObs[i]), sigmaPK);

  for (i in 1:nObsPD)
    pdObs_pred[i] = lognormal_rng(log(pdHatObs[i]), sigmaPD);
}
