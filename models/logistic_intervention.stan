
data {
  int<lower=0> Nsubj;  // number of subjects
  int<lower=0> Nscen;  // number of cases
  int<lower=0> N;  // number of observations
  int<lower=0> P;  // number of fixed + random effect regressors
  matrix[N, P] X;  // design matrix for fixed effects
  int<lower=0> Scen[N];  // subject corresponding to each rating
  int<lower=0> Subj[N];  // case corresponding to each rating
  int Y[N]; // guilt judgement
  int Z[N]; // intervention dummy
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
  
  vector[2] mu_lambda; //intervention effect
  vector<lower=0>[2] sigma_lambda;
  cholesky_factor_corr[2] L_lambda;
  vector[2] lambda_subj_raw[Nsubj];
  
}

transformed parameters {
  vector[P] beta_scen[Nscen];  // scenario effects
  vector[P] beta_subj[Nsubj];  // individual effects
  vector[Nscen] alpha_scen = sigma_alpha_scen * alpha_scen_raw;
  vector[Nsubj] alpha_subj = sigma_alpha_subj * alpha_subj_raw;
  vector[2] lambda_subj[Nsubj]; // individual intervention effects
  real log_lik[N];
  
  //random effects
  for (i in 1:Nscen)
  beta_scen[i] = sigma_beta_scen .* beta_scen_raw[i];
  for (i in 1:Nsubj) {
    beta_subj[i] = sigma_beta_subj .* beta_subj_raw[i];
    lambda_subj[i] = mu_lambda + diag_pre_multiply(sigma_lambda,L_lambda) * lambda_subj_raw[i];

  }
  //linear predictor  
  for (i in 1:N) {
    real alpha = mu_alpha + alpha_subj[Subj[i]] + alpha_scen[Scen[i]];
    vector[P] beta = mu_beta + beta_subj[Subj[i]] + beta_scen[Scen[i]];
    real eta = alpha + Z[i]*lambda_subj[Subj[i]][1] + X[i]*beta + sum(X[i])*Z[i]*lambda_subj[Subj[i]][2];
    log_lik[i] = bernoulli_logit_lpmf(Y[i] | eta);
  }
}

model {
  
  target += sum(log_lik);
  
  mu_alpha ~ normal(0, 2.5);
  sigma_alpha_scen ~ normal(0, 2);
  sigma_alpha_subj ~ normal(0, 2);
  
  mu_beta ~ normal(0, 2.5);
  sigma_beta_scen ~ normal(0, 2);
  sigma_beta_subj ~ normal(0, 2);
  
  mu_lambda ~ normal(0,2.5);
  sigma_lambda ~ normal(0,2);
  
  alpha_subj_raw ~ normal(0,1);
  alpha_scen_raw ~ normal(0,1);
  
  for (i in 1:Nsubj) {
    beta_subj_raw[i] ~ normal(0, 1);
    lambda_subj_raw[i] ~ normal(0,1);
  }
  
  for (i in 1:Nscen) beta_scen_raw[i] ~ normal(0, 1);
}

generated quantities {
    real rho_lambda = (L_lambda * L_lambda')[2,1];
}
