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
  int<lower=0> Nsubj;  // number of subjects
  int<lower=0> Nscen;  // number of cases
  int<lower=0> N;  // number of observations
  int<lower=0> P;  // number of fixed + random effect regressors
  real<lower=L, upper=U> Y[N];  // ratings
  matrix[N, P] X;  // design matrix for fixed + random effects
  int Z[N]; //intervention
  real W[Nsubj]; //intervention interaction
  int<lower=0> Subj[N];  // subject corresponding to each rating
  int<lower=0> Scen[N];  // case corresponding to each rating
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
  
  vector[2] mu_lambda; //intervention effect
  vector<lower=0>[2] sigma_lambda;
  cholesky_factor_corr[2] L_lambda;
  vector[2] lambda_subj_raw[Nsubj];
  
  vector[2] gamma;
  
  real mu_scale;
  real<lower=0> sigma_scale;
  vector[Nsubj] scale_subj_raw; //scale of internal transfer function
}

transformed parameters {
  vector[P] beta_scen[Nscen];  // scenario effects
  vector[P] beta_subj[Nsubj];  // individual effects
  vector[2] lambda_subj[Nsubj]; // individual intervention effects
  vector[Nscen] alpha_scen = sigma_alpha_scen * alpha_scen_raw;
  vector[Nsubj] alpha_subj = sigma_alpha_subj * alpha_subj_raw;
  vector[Nsubj] scale_subj = exp(mu_scale + sigma_scale * scale_subj_raw);

  // random effects
  for (i in 1:Nscen)
    beta_scen[i] = sigma_beta_scen .* beta_scen_raw[i];
  for (i in 1:Nsubj) {
    beta_subj[i] = sigma_beta_subj .* beta_subj_raw[i];
    lambda_subj[i] = diag_pre_multiply(sigma_lambda,L_lambda) * lambda_subj_raw[i];
  }
  
}

model {
  
  for (i in 1:N) {
    
    real alpha = mu_alpha + alpha_subj[Subj[i]] + alpha_scen[Scen[i]];
    vector[P] beta = mu_beta + beta_subj[Subj[i]] + beta_scen[Scen[i]];
    vector[2] lambda = mu_lambda + lambda_subj[Subj[i]] + gamma*W[Subj[i]];
    real eta = alpha + Z[i]*lambda[1] + X[i]*beta + sum(X[i])*Z[i]*lambda[2];

    if (Y[i] == L)
      target += lcdf(inv_cdf(I)*scale_subj[Subj[i]], eta, 1);
    else if (Y[i] == U)
      target += lccdf(inv_cdf(1-I)*scale_subj[Subj[i]], eta, 1);
    else
      target += log_diff_exp(lcdf(inv_cdf(Q[i]+I)*scale_subj[Subj[i]], eta, 1), lcdf(inv_cdf(Q[i]-I)*scale_subj[Subj[i]], eta, 1));
    }
  
  mu_alpha ~ normal(0,2.5);
  sigma_alpha_scen ~ normal(0,2);
  sigma_alpha_subj ~ normal(0,2);
  
  mu_beta ~ normal(0,2.5);
  sigma_beta_scen ~ normal(0,2);
  sigma_beta_subj ~ normal(0,2);
  
  mu_lambda ~ normal(0,2.5);
  sigma_lambda ~ normal(0,2);
  
  gamma ~ normal(0,1);
  
  alpha_subj_raw ~ normal(0,1);
  alpha_scen_raw ~ normal(0,1);
  
  for (i in 1:Nsubj) {
    beta_subj_raw[i] ~ normal(0,1);
    lambda_subj_raw[i] ~ normal(0,1);
  }
  for (i in 1:Nscen)
    beta_scen_raw[i] ~ normal(0,1);
  
  mu_scale ~ normal(0,1);
  sigma_scale ~ normal(0,1);
  scale_subj_raw ~ normal(0,1);
}

generated quantities {
  real mu_alpha_resp = Phi(mu_alpha/exp(mu_scale + pow(sigma_scale,2)/2))*(U-L) + L;
  vector[P] mu_beta_resp;
  real mu_lambda_alpha = (Phi(mu_lambda[1]/exp(mu_scale + pow(sigma_scale,2)/2)) - 0.5)*(U-L) + L;
  real mu_lambda_beta = (Phi(mu_lambda[2]/exp(mu_scale + pow(sigma_scale,2)/2)) - 0.5)*(U-L) + L;
  real rho_lambda = (L_lambda * L_lambda')[2,1];

  for (p in 1:P) mu_beta_resp[p] = (Phi(mu_beta[p]/exp(mu_scale + pow(sigma_scale,2)/2)) - 0.5)*(U-L) + L;
}
