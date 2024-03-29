functions {
  int sign(real x) {
    return x < 0 ? -1 : x > 0;
  }

  real[] etaize(matrix X, int[] Subj, int[] Scen, real mu_alpha, 
                  vector alpha_subj, vector alpha_scen, vector mu_beta, vector[] beta_subj, vector[] beta_scen, real gamma) {
    
    int N = num_elements(Subj);
    real eta[N];
    for (i in 1:num_elements(Subj)) {
        eta[i] = mu_alpha + alpha_scen[Scen[i]] + alpha_subj[Subj[i]] + X[i]*(mu_beta + beta_scen[Scen[i]] + beta_subj[Subj[i]]);
        eta[i] = sign(eta[i]) * pow(abs(eta[i]),gamma);
    }
    return eta;
  }
  
  real inv_cdf(real p) {
    return inv_Phi(p);
  }
  
  real lcdf(real q, real mu, real sigma) {
    return logistic_lcdf(q | mu, sigma);
  }
  
  real lccdf(real q, real mu, real sigma) {
    return logistic_lccdf(q | mu, sigma);
  }
  
  real[] Yhatify(real[] eta, vector scale, int[] Subj, int L, int U, int D) {
    int N = num_elements(eta);
    real Yhat[N];
    int M = (U-L)/D - 1;
    real I = D/(2.*(U-L));

    for (i in 1:N) {
      Yhat[i] = exp(lccdf(inv_cdf(1-I)*scale[Subj[i]], eta[i], 1));
      for (j in 1:M)
        Yhat[i] = Yhat[i] + I*2*j*exp(log_diff_exp(lcdf(inv_cdf(I*(2*j+1))*scale[Subj[i]], eta[i], 1),
                                  lcdf(inv_cdf(I*(2*j-1))*scale[Subj[i]], eta[i], 1)));
      Yhat[i] = Yhat[i]*(U-L) + L;
    }
    return Yhat;
  }
}

data {
  int L;  // lower censoring
  int U;  // upper censoring
  int D; // distance between consecutive points on scale
  int<lower=0> Nsubj;  // number of subjects
  int<lower=0> Nscen;  // number of cases
  int<lower=0> N;  // number of observations
  int<lower=0> P;  // number of fixed + random effect regressors
  real<lower=L, upper=U> Y[N];  // ratings
  matrix[N, P] X;  // design matrix for fixed + random effects
  int<lower=0> Subj[N];  // subject corresponding to each rating
  int<lower=0> Scen[N];  // case corresponding to each rating
  // int<lower=1> K;    //components of mixture transfer function
}

transformed data {
  real I = D/(2.*(U-L));
  real<lower=0,upper=1> Q[N];

  for (i in 1:N)
    Q[i] = (Y[i]-L)/(U-L);
}

parameters {
  real mu_alpha;
  real<lower=0> sigma_alpha_scen;
  real<lower=0> sigma_alpha_subj;
  
  vector[Nscen] alpha_scen_raw;
  vector[Nsubj] alpha_subj_raw;

  vector[P] mu_beta;
  vector<lower=0>[P] sigma_beta_scen;
  vector<lower=0>[P] sigma_beta_subj;

  vector[P] beta_scen_raw[Nscen];  // scenario effects
  vector[P] beta_subj_raw[Nsubj];  // subject residual effects
  
  real mu_scale;
  real<lower=0> sigma_scale;
  vector[Nsubj] scale_subj_raw; //scale of internal transfer function
  
  real<lower=0,upper=2> gamma; //lol
}

transformed parameters {
  vector[P] beta_scen[Nscen];  // scenario effects
  vector[P] beta_subj[Nsubj];  // individual effects
  vector[Nscen] alpha_scen = sigma_alpha_scen * alpha_scen_raw;
  vector[Nsubj] alpha_subj = sigma_alpha_subj * alpha_subj_raw;
  vector[Nsubj] scale_subj = exp(mu_scale + sigma_scale * scale_subj_raw);
  real eta[N];
  real log_lik[N];

  //random effects
  for (i in 1:Nscen) 
    beta_scen[i] = sigma_beta_scen .* beta_scen_raw[i];
  for (i in 1:Nsubj)
    beta_subj[i] = sigma_beta_subj .* beta_subj_raw[i];
  
  // get linear predictor
  eta = etaize(X, Subj, Scen, mu_alpha, alpha_subj, alpha_scen, mu_beta, beta_subj, beta_scen, gamma);

  // eval log lik
  for (i in 1:N) {
    if (Y[i] == L)
      log_lik[i] = lcdf(inv_cdf(I)*scale_subj[Subj[i]], eta[i], 1);
    else if (Y[i] == U)
      log_lik[i] = lccdf(inv_cdf(1-I)*scale_subj[Subj[i]], eta[i], 1);
    else
      log_lik[i] = log_diff_exp(lcdf(inv_cdf(Q[i]+I)*scale_subj[Subj[i]], eta[i], 1), lcdf(inv_cdf(Q[i]-I)*scale_subj[Subj[i]], eta[i], 1));
    }
}

model {
  
  for (i in 1:N) target += log_lik[i];
  
  mu_alpha ~ normal(0,2.5);
  sigma_alpha_scen ~ normal(0,2);
  sigma_alpha_subj ~ normal(0,2);
  
  mu_beta ~ normal(0,2.5);
  sigma_beta_scen ~ normal(0,2);
  sigma_beta_subj ~ normal(0,2);
  
  alpha_subj_raw ~ normal(0,1);
  alpha_scen_raw ~ normal(0,1);
  
  for (i in 1:Nsubj)
    beta_subj_raw[i] ~ normal(0,1);
  for (i in 1:Nscen)
    beta_scen_raw[i] ~ normal(0,1);
  
  mu_scale ~ normal(0,1);
  sigma_scale ~ normal(0,1);
  scale_subj_raw ~ normal(0,1);
  gamma ~ normal(1,0.1);
}

generated quantities {
  real Yhat[N] = Yhatify(eta, scale_subj, Subj, L, U, D);
  real mu_alpha_resp = Phi(mu_alpha/exp(mu_scale + pow(sigma_scale,2)/2))*(U-L) + L;
  vector[P] mu_beta_resp;

  for (p in 1:P) mu_beta_resp[p] = (Phi(mu_beta[p]/exp(mu_scale + pow(sigma_scale,2)/2)) - 0.5)*(U-L) + L;
}
