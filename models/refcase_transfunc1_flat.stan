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
  int<lower=0> NX[N];
}

transformed data {
  real I;
  real<lower=0,upper=1> Q[N];
  int<lower=-1,upper=1> cens[N];
  real QU;
  real QL;
  int Ncond = max(Cond);

  
  for (i in 1:N) {
    if (Y[i] > (U-1.*D)) {
      cens[i] = 1;
    } else if (Y[i] < (L+1.*D)) {
      cens[i] = -1;
    } else {
      cens[i] = 0;
    }
  }
  
  I = D/(2.*(U-L));
  QU = (U-1.*D+D-L)/(U-L);
  QL = (1.*D-D)/(U-L);
  for (i in 1:N)
    Q[i] = (Y[i]-L)/(U-L);
}

parameters {
  real alpha;
  vector[P] beta;
  real<lower=0> scale; //scale of internal transfer function
  real<lower=0,upper=1> tau;
  real<lower=0> nu[Ncond];
  real<lower=0> xi[Ncond];
  real W0;
}

transformed parameters {
  real eta[N];
  real Wbar[N];
  
  for (i in 1:N) {
    real delta;
    // real W0 = 0;
    real j;
    
    if (Trial[i] == 1) {
      Wbar[i] = W0;
    } else {
      if (NX[i-1] == 0) delta = -Wbar[i-1];
      else delta = (X[i-1]*beta)./NX[i-1] - Wbar[i-1];
      Wbar[i] = Wbar[i-1] + tau*delta;
    }
    eta[i] = alpha + X[i]*(beta-Wbar[i]) + (nu[Cond[i]] + xi[Cond[i]]*NX[i])*Wbar[i];
  }
  
}
model {
  for (i in 1:N) {
    if (cens[i] == 0)
      target += log_diff_exp(logistic_lcdf(inv_Phi(Q[i]+I)*scale | eta[i], 1), logistic_lcdf(inv_Phi(Q[i]-I)*scale | eta[i], 1));
    else if (cens[i] == -1)
      target += logistic_lcdf(inv_Phi(QL+I)*scale | eta[i], 1);
    else if (cens[i] == 1)
      target += logistic_lccdf(inv_Phi(QU-I)*scale | eta[i], 1);
  }
}

generated quantities {
  real baseline[N];
  real exceff[N];
  real inceff[N];
  real ambeff[N];
  real combined[N];
  
  for (i in 1:N) {
    baseline[i] = alpha + nu[Cond[i]]*Wbar[i];
    exceff[i] = (beta[1] + beta[4] + beta[7] + beta[10])/4 + (xi[Cond[i]]-1)*Wbar[i];
    ambeff[i] = (beta[2] + beta[5] + beta[8])/3 + (xi[Cond[i]]-1)*Wbar[i];
    inceff[i] = (beta[3] + beta[6] + beta[9] + beta[11])/4 + (xi[Cond[i]]-1)*Wbar[i];
    combined[i] = mean(beta) + (xi[Cond[i]]-1)*Wbar[i];
  }
}
