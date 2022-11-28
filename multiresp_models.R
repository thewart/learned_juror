source("~/code/casereveal/analysis/process_data.R")
source("~/code/casereveal/analysis/prep_stan.R")
ldat <- makelegaldat(scendat,subjdat,clickdat)
standat <- makestandat(ldat,binresp=T)
# standat$Y <- ldat[,.(as.numeric(subjguilt),as.numeric(legalguilt))]
# standat$D <- 2
standat$Y <- ldat[,.(as.numeric(subjguilt),as.numeric(likelyguilt),as.numeric(legalguilt))]
standat$D <- 3


# standat$Y <- ldat[,.(legalguilt)]
# standat <- c(standat,D=1,D2=2)

# standat$X <- standat$X[,-1]
# standat$P <- standat$P - 1
standat$D2 <- 2^standat$D
qeb <- stan_model("~/code/legalmodels/qeb_lapse.stan")
fit <- sampling(qeb,standat,chains=4,iter=400,warmup=200,
                pars=c("mu_alpha","mu_beta","mu_f2","mu_eps","sigma_alpha_subj","sigma_alpha_scen",
                       "sigma_beta_subj","sigma_beta_scen","sigma_f2_subj","sigma_f2_scen","sigma_eps","log_lik"))

standat$Y <- t(standat$Y)
qeb <- stan_model("~/code/legalmodels/qeb_lapse_allreg.stan");
fit <- sampling(qeb,standat,chains=4,iter=400,warmup=200,
                pars=c("mu_alpha","mu_beta","mu_eps","sigma_alpha",
                       "sigma_beta","sigma_eps","Rho_alpha","Rho_beta"))


m0 <- stan_model("~/code/legalmodels/logistic_lapse.stan")
f0 <- sampling(m0,st0,chains=1,iter=300,warmup=200,
                pars=c("beta_mu","sigma_subj","sigma_scen"))


# standat$X <- standat$X[,-1]
# standat$P <- standat$P - 1
YO <- ldat[,subjguilt+legalguilt+1]
YO[ldat[,subjguilt==F & legalguilt==T]] <- 0
standat$Y <- YO
standat$ps <- 1e6
ord <- stan_model("~/code/legalmodels/ordreg_lapse.stan")
fot <- sampling(ord,standat,chains=4,iter=400,warmup=150)
                # pars=c("mu_alpha","sigma_alpha","mu_theta","sigma_theta","mu_legalgap","sigma_legalgap","mu_lapserate","sigma_lapserate",
                #        "rho_alpha_theta","mu_beta","mu_eps","sigma_eps","sigma_subj","sigma_scen","alpha","theta","log_lik"))
