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
  int<lower=0> m[N];
}

transformed data {
  int Ncond = max(Cond);
  real I = D/(2.*(U-L));
  real<lower=0,upper=1> Q[N];

  for (i in 1:N)
    Q[i] = (Y[i]-L)/(U-L);
}

parameters {
  real alpha;
  vector[P] beta;
  real<lower=0> scale; //scale of internal transfer function
  real<lower=0,upper=1> tau;
  real<lower=0> nu_0[Ncond];
  real<lower=0> nu_1[Ncond];
  real W0;
}

transformed parameters {
  real eta[N];
  real Wbar[N];
  
  for (i in 1:N) {
    real nhat;
    real delta;
    // real W0 = 0;
    real j;
    
    if (Trial[i] == 1) {
      Wbar[i] = W0;
    } else {
      if (m[i-1] == 0) delta = -Wbar[i-1];
      else delta = (X[i-1]*beta)./m[i-1] - Wbar[i-1];
      Wbar[i] = Wbar[i-1] + tau*delta;
    }
    nhat = nu_0[Cond[i]] + nu_1[Cond[i]]*m[i]
    eta[i] = alpha + X[i]*beta + (nhat-m[i])*Wbar[i];
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
    baseline[i] = alpha + nu_0[Cond[i]]*Wbar[i];
    exceff[i] = (beta[1] + beta[4] + beta[7] + beta[10])/4 + (nu_1[Cond[i]]-1)*Wbar[i];
    ambeff[i] = (beta[2] + beta[5] + beta[8])/3 + (nu_1[Cond[i]]-1)*Wbar[i];
    inceff[i] = (beta[3] + beta[6] + beta[9] + beta[11])/4 + (nu_1[Cond[i]]-1)*Wbar[i];
    combined[i] = mean(beta) + (nu_1[Cond[i]]-1)*Wbar[i];
  }
}
