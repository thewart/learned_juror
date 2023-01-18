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
  // int<lower=0> Nscen;  // number of cases
  int<lower=0> N;  // number of observations
  int<lower=0> P;  // number of fixed + random effect regressors
  real<lower=L, upper=U> Y[N];  // ratings
  matrix[N, P] X;  // design matrix for fixed + random effects
  vector[N] Z; //intervention
  int<lower=0> Subj[N];  // subject corresponding to each rating
  // int<lower=0> Scen[N];  // case corresponding to each rating
}

transformed data {
  real I = D/(2.*(U-L));
  real<lower=0,upper=1> Q[N];

  for (i in 1:N)
    Q[i] = (Y[i]-L)/(U-L);
}

parameters {
  real mu_alpha;
  // real<lower=0> sigma_alpha_scen;
  real<lower=0> sigma_alpha_subj;
  
  // vector[Nscen] alpha_scen_raw;
  vector[Nsubj] alpha_subj_raw;

  vector[P] mu_beta;
  // vector<lower=0>[P] sigma_beta_scen;
  vector<lower=0>[P] sigma_beta_subj;

  // vector[P] beta_scen_raw[Nscen];  // scenario effects
  vector[P] beta_subj_raw[Nsubj];  // subject residual effects
  
  vector[2] mu_lambda; //intervention effect
  vector<lower=0>[2] sigma_lambda;
  corr_matrix[2] Omega_lambda;
  vector[2] lambda_subj_raw[[Nsubj]];
  
  real mu_scale;
  real<lower=0> sigma_scale;
  vector[Nsubj] scale_subj_raw; //scale of internal transfer function
}

transformed parameters {
  // vector[P] beta_scen[Nscen];  // scenario effects
  vector[P] beta_subj[Nsubj];  // individual effects
  vector[2] lambda[Nsubj]
  // vector[Nscen] alpha_scen = sigma_alpha_scen * alpha_scen_raw;
  vector[Nsubj] alpha_subj = sigma_alpha_subj * alpha_subj_raw;
  vector[Nsubj] scale_subj = exp(mu_scale + sigma_scale * scale_subj_raw);

  //random effects
  // for (i in 1:Nscen) 
    // beta_scen[i] = sigma_beta_scen .* beta_scen_raw[i];
  for (i in 1:Nsubj) {
    beta_subj[i] = sigma_beta_subj .* beta_subj_raw[i];
    lambda_subj[i] = quad_form_diag(Omega_lambda,sigma_lambda) * lambda_subj_raw[i];
  }
  
}

model {

  for (i in 1:N) {
    real alpha = mu_alpha + alpha_subj[Subj[i]];
    real beta = mu_beta + beta_subj[Subj[i]];
    real lambda = mu_lambda + lambda_subj[i];
    real eta = alpha + Z[i]*lambda[1] + X[i]*beta + X[i]*Z[i]*lambda[2];
    
    if (Y[i] == L)
      target += lcdf(inv_cdf(I)*scale_subj[Subj[i]], eta, 1);
    else if (Y[i] == U)
      target += lccdf(inv_cdf(1-I)*scale_subj[Subj[i]], eta, 1);
    else
      target += log_diff_exp(lcdf(inv_cdf(Q[i]+I)*scale_subj[Subj[i]], eta, 1), lcdf(inv_cdf(Q[i]-I)*scale_subj[Subj[i]], eta, 1));
  }
    
  mu_alpha ~ normal(0,2.5);
  // sigma_alpha_scen ~ normal(0,2);
  sigma_alpha_subj ~ normal(0,2);
  
  mu_beta ~ normal(0,2.5);
  // sigma_beta_scen ~ normal(0,2);
  sigma_beta_subj ~ normal(0,2);
  
  // alpha_scen_raw ~ normal(0,1);
  alpha_subj_raw ~ normal(0,1);
  
  mu_lambda ~ normal(0,2.5);
  sigma_lambda_subj ~ normal(2);
  
  for (i in 1:Nsubj) {
    beta_subj_raw[i] ~ normal(0,1);
    lambda_subj_raw[i] ~ normal(0,1);
  }
  // for (i in 1:Nscen)
    // beta_scen_raw[i] ~ normal(0,1);
  
  mu_scale ~ normal(0,1);
  sigma_scale ~ normal(0,1);
  scale_subj_raw ~ normal(0,1);
}

