functions {
  real[] etaize(matrix X, int[] Subj, int[] Scen, int[] Cond, real[] mu_alpha, 
  vector alpha_subj, vector alpha_scen, vector[] mu_beta, vector[] beta_subj, vector[] beta_scen) {
    
    int N = num_elements(Subj);
    real eta[N];
    for (i in 1:N) {
      eta[i] = mu_alpha[Cond[i]] + alpha_scen[Scen[i]] + alpha_subj[Subj[i]] + X[i]*(mu_beta[Cond[i]] + beta_scen[Scen[i]] + beta_subj[Subj[i]]);
    }
    return eta;
  }
}

data {
  int<lower=0> Nsubj;  // number of subjects
  int<lower=0> Nscen;  // number of cases
  int<lower=0> Ncond;  // number of conditions
  int<lower=0> N;  // number of observations
  int<lower=0> P;  // number of fixed + random effect regressors
  matrix[N, P] X;  // design matrix for fixed effects
  int<lower=0> Scen[N];  // subject corresponding to each response
  int<lower=0> Subj[N];  // case corresponding to each response
  int<lower=0> Cond[N];  // condition for each response
  int Y[N]; // guilt judgement
}

parameters {
  
  real mu_alpha[Ncond];
  real<lower=0> sigma_alpha_scen;
  real<lower=0> sigma_alpha_subj;
  
  vector[Nscen] alpha_scen_raw;
  vector[Nsubj] alpha_subj_raw;
  
  vector[P] mu_beta[Ncond];
  vector<lower=0>[P] sigma_beta_scen;
  vector<lower=0>[P] sigma_beta_subj;
  
  vector[P] beta_scen_raw[Nscen];  // scenario effects
  vector[P] beta_subj_raw[Nsubj];  // subject residual effects
  
}

transformed parameters {
  vector[P] beta_scen[Nscen];  // scenario effects
  vector[P] beta_subj[Nsubj];  // individual effects
  vector[Nscen] alpha_scen = sigma_alpha_scen * alpha_scen_raw;
  vector[Nsubj] alpha_subj = sigma_alpha_subj * alpha_subj_raw;
  real log_lik[N];
  
  //random effects
  for (i in 1:Nscen) beta_scen[i] = sigma_beta_scen .* beta_scen_raw[i];
  for (i in 1:Nsubj) beta_subj[i] = sigma_beta_subj .* beta_subj_raw[i];
  
  //linear predictor  
  {
    real eta[N] = etaize(X, Subj, Scen, Cond, mu_alpha, alpha_subj, alpha_scen, mu_beta, beta_subj, beta_scen);
    for (i in 1:N) log_lik[i] = bernoulli_logit_lpmf(Y[i] | eta[i]);
  }
}

model {
  
  
  target += sum(log_lik);
  
  for (c in 1:Ncond) {
    mu_alpha[c] ~ normal(0,2.5);
    mu_beta[c] ~ normal(0,2.5);
  }
  
  sigma_alpha_scen ~ normal(0, 2);
  sigma_alpha_subj ~ normal(0, 2);
  
  sigma_beta_scen ~ normal(0, 2);
  sigma_beta_subj ~ normal(0, 2);
  
  alpha_subj_raw ~ normal(0,1);
  alpha_scen_raw ~ normal(0,1);
  for (i in 1:Nsubj) beta_subj_raw[i] ~ normal(0, 1);
  for (i in 1:Nscen) beta_scen_raw[i] ~ normal(0, 1);
}
