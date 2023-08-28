functions {
  real inv_cdf(real p) {
    return inv_Phi(p);
  }
  
  real lcdf(real q, real mu, real sigma) {
    return logistic_lcdf(q | mu, sigma);
  }
  
  real lccdf(real q, real mu, real sigma) {
    return logistic_lccdf(q | mu, sigma);
  }
}

data {
  int L;  // lower censoring
  int U;  // upper censoring
  int D; // distance between consecutive points on scale
  int<lower=0> N;  // number of observations
  int<lower=0> P;  // number of fixed + random effect regressors
  real<lower=L,upper=U> Y[N];  // ratings
  matrix[N,P] X;  // design matrix for fixed + random effects
  int<lower=0> Trial[N];
  int<lower=1> Cond[N];
  int<lower=0> m[N]; // amount of evidence on each trial
}

transformed data {
  int Ncond = max(Cond);
  real I = D/(2.*(U-L));
  real<lower=0,upper=1> Q[N];
  real N_exc[N];
  real N_inc[N];
  real N_amb[N];
  
  for (i in 1:N) {
    Q[i] = (Y[i]-L)/(U-L);
    N_exc[i] = X[i][1] + X[i][4] + X[i][7] + X[i][10];
    N_inc[i] = X[i][3] + X[i][6] + X[i][9] + X[i][11];
    N_amb[i] = X[i][2] + X[i][5] + X[i][8];
  }
}

parameters {
  // real alpha;
  vector[P] beta;
  real<lower=0> scale; //scale of internal transfer function
  real<lower=0,upper=1> tau;
  real<lower=0> nu_0; //nhat intercept
  real<lower=0> nu_1; //nhat slope
  real psi_0; //truncation prob intercept
  real psi_inc; //truncation prob slope
  real psi_amb; //truncation prob slope
  real psi_exc; //truncation prob slope
  real gamma; //truncation bonus
}

transformed parameters {
  real eta[N];
  real Wbar[N];
  real trunc_bonus[N];
  
  {
    real N_exc_cum;
    real N_inc_cum;
    real N_amb_cum;
    for (i in 1:N) {
      real W0 = mean(beta);
      real delta;
      real nhat;
      real j;
      
      if (Trial[i] == 1) {
        Wbar[i] = W0;
        N_exc_cum = 0;
        N_inc_cum = 0;
        N_amb_cum = 0;
      } else {
        if (m[i-1] == 0) delta = -Wbar[i-1];
        else delta = (X[i-1]*beta)./m[i-1] - Wbar[i-1];
        Wbar[i] = Wbar[i-1] + tau*delta;
        N_exc_cum = N_exc_cum + N_exc[i-1];
        N_inc_cum = N_inc_cum + N_inc[i-1];
        N_amb_cum = N_amb_cum + N_amb[i-1];
      }
      
      trunc_bonus[i] = inv_logit(psi_0 + psi_inc*N_inc_cum + psi_exc*N_exc_cum) * gamma;
      nhat = nu_0 + nu_1*m[i];
      eta[i] = X[i]*beta + (nhat-m[i])*(Wbar[i] + trunc_bonus[i]);
    }
  }
}

model {
  for (i in 1:N) {
    if (Y[i] == L)
    target += lcdf(inv_cdf(I)*scale, eta[i], 1);
    else if (Y[i] == U)
    target += lccdf(inv_cdf(1-I)*scale, eta[i], 1);
    else
    target += log_diff_exp(lcdf(inv_cdf(Q[i]+I)*scale, eta[i], 1), lcdf(inv_cdf(Q[i]-I)*scale, eta[i], 1));
  }
}

generated quantities {
  real baseline[N];
  real exceff[N];
  real inceff[N];
  real ambeff[N];
  real combined[N];
  
  for (i in 1:N) {
    baseline[i] = nu_0*(Wbar[i] + trunc_bonus[i]);
    exceff[i] = (beta[1] + beta[4] + beta[7] + beta[10])/4 + (nu_1-1)*(Wbar[i] + trunc_bonus[i]);
    ambeff[i] = (beta[2] + beta[5] + beta[8])/3 + (nu_1-1)*(Wbar[i] + trunc_bonus[i]);
    inceff[i] = (beta[3] + beta[6] + beta[9] + beta[11])/4 + (nu_1-1)*(Wbar[i] + trunc_bonus[i]);
    combined[i] = mean(beta) + (nu_1-1)*(Wbar[i] + trunc_bonus[i]);
  }
}
