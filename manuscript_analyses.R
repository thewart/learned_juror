model <- stan_model("~/code/legalmodels/transfunc1.stan")
intmodel <- stan_model("~/code/legalmodels/transfunc1_interactions.stan")
bardmodel <- stan_model("~/code/legalmodels/logistic.stan")
basepars <- c("mu_alpha","mu_beta","sigma_alpha_subj","sigma_alpha_scen","sigma_beta_subj","sigma_beta_scen",
              "alpha_scen","beta_scen","alpha_subj","beta_subj","log_lik")
rpars <- c("mu_alpha_resp","mu_beta_resp","mu_scale","sigma_scale","scale_subj","Yhat")
intpars <- c("sigma_mu_lambda","sigma_lambda_subj","mu_lambda","lambda_subj")
evcond <- c("balanced","credible","defenseless")
ratecond <- c("without","with")

source("~/code/casereveal/analysis/miscfunctions.R")
source("~/code/casereveal/analysis/prep_stan.R")

#### Experiment 1 #####

whichtask <- "exculpatory_burdenofproof"
source("~/code/casereveal/analysis/process_data.R")
ldat <- makelegaldat(scendat,subjdat,clickdat)

##### split fits #####

ldat_2019 <- ldat[uid %in% subjdat[,as.POSIXct(start) %>% format(format="%Y"),by=uid][V1=="2019",uid]]
ldat_2020 <- ldat[uid %in% subjdat[,as.POSIXct(start) %>% format(format="%Y"),by=uid][V1=="2020",uid]]
standat_2019 <- makestandat(ldat_2019, interact = T)
standat_2020 <- makestandat(ldat_2020, interact = T)

rfit_2020 <- sampling(model,standat_2020,chains=4,init=transfunc_init,iter=450,warmup=200,
                      pars=c(basepars,rpars,"Yhat","log_lik"))

rfit_2019 <- sampling(model,standat_2019,chains=4,init=transfunc_init,iter=450,warmup=200,
                      pars=c(basepars,rpars,"Yhat","log_lik"))

rfitint_2019 <- sampling(intmodel,standat_2019,chains=4,init=transfunc_init,iter=450,warmup=200,
                        pars=c(basepars,rpars,intpars,"Yhat","log_lik"))

rfitint_2020 <- sampling(intmodel,standat_2020,chains=4,init=transfunc_init,iter=450,warmup=200,
                         pars=c(basepars,rpars,intpars,"Yhat","log_lik"))

##### joint fits #####

bfit <- sampling(bardmodel,makestandat(ldat,binresp=T),chains=4,iter=1000,warmup=500,pars=basepars)

standat <- makestandat(ldat, interact = T)
rfit <- sampling(model,standat,chains=4,init=transfunc_init,iter=1000,warmup=500,
                      pars=c(basepars,rpars))
rfitint <- sampling(intmodel,standat,chains=4,init=transfunc_init,iter=1000,warmup=500,
                        pars=c(basepars,rpars,intpars))


##### PPC -- Case strength rating distribution ####
expose_stan_functions(rfit)
ob <- c(standat,extract(rfit))
yrep <- postpredsample(ob)

##### Interactions ####
looout <- loo_compare(loo(rfit),loo(rfitint))
expose_stan_functions(rfitint)
cyhat <- postpredeta_cfact(c(extract(rfitint),standat),c("mu_lambda","lambda_subj"),
                              binint=T,list_out=T) %>% postpredyhat()

##### Scenario weights ####
expose_stan_functions(rfit)
beta_scen_resp <- get_beta_scen_resp(rfit)
alpha_scen_resp <- get_alpha_scen_resp(rfit)

##### Marginal impact of evidence ####
library(rstanarm)
expose_stan_functions(rfit)
beta <- extract(rfit,"mu_beta")[[1]]
alpha <- extract(rfit,"mu_alpha")[[1]]
scale <- extract(rfit,"mu_scale")[[1]]
scalescale <- extract(rfit,"sigma_scale")[[1]]

bdat <- merge(makebalancedat(ldat),ldat,by=c("uid","scenario"))
bdat[,evbal := n_inculp - n_exculp]
foo <- bdat[,.(n_exculp,n_inculp,evbal)]
juh <- bdat[,.(rating,uid,evbal)]
marginaldiff <- data.table()

standat <- makestandat(ldat)

for (i in -3:3) {
  fooi <- unique(foo[evbal==i])
  inde <- vector(length=nrow(juh))
  indi <- vector(length=nrow(juh))
  
  for (j in 1:nrow(fooi)) {
    inde[foo[,(n_exculp == fooi[j,n_exculp+1]) & (n_inculp == fooi[j,n_inculp])]] <- T
    indi[foo[,n_inculp == fooi[j,n_inculp+1] & n_exculp == fooi[j,n_exculp]]] <- T
  }
  
  juh$plus_inc <- indi
  juh$plus_exc <- inde
  
  dat_inc <- juh[evbal==i | plus_inc==T]
  guh <- stan_lmer(scale(rating) ~ 1 + plus_inc + (1 + plus_inc | uid), data=dat_inc,iter=400)$stanfit
  # control=lmerControl(check.nobs.vs.nRE = "ignore"))
  # se <- sqrt(vcov(guh)[2,2])
  post_inc <- dat_inc[,extract(guh,"plus_incTRUE")[[1]] * sd(rating-mean(rating))] %>% post_summary_dt()
  marginaldiff <- rbind(marginaldiff,cbind(evbal=i,post_inc,type="Inculpatory",source="Behavior"))
  
  dat_exc <- juh[evbal==i | plus_exc==T]
  guh <- stan_lmer(scale(rating) ~ 1 +plus_exc + (1 + plus_exc | uid), data=dat_exc, iter=400)$stanfit
  # control=lmerControl(check.nobs.vs.nRE = "ignore"))
  # se <- sqrt(vcov(guh)[2,2])
  post_exc <- dat_exc[,extract(guh,"plus_excTRUE")[[1]] * sd(rating-mean(rating))] %>% post_summary_dt()
  marginaldiff <- rbind(marginaldiff,cbind(evbal=i,post_exc,type="Exculpatory",source="Behavior"))
  
  evbal_configs <- unique(standat$X[foo[,evbal==i],])
  eta0 <- matrix(alpha,nrow=2000,ncol=nrow(evbal_configs)) + beta %*% t(evbal_configs)
  yhat0 <- array(dim = dim(eta0))
  for (k in 1:nrow(evbal_configs))
    yhat0[,k] <- Yhatify(eta0[,k],exp(scale + scalescale^2/2),1:2000,standat$L,standat$U,standat$D)
  yhat0 <- rowMeans(yhat0)
  
  plusexc_configs <- unique(standat$X[inde,])
  etae <- matrix(alpha,nrow=2000,ncol=nrow(plusexc_configs)) + beta %*% t(plusexc_configs)
  yhate <- array(dim = dim(etae))
  for (k in 1:nrow(plusexc_configs))
    yhate[,k] <- Yhatify(etae[,k],exp(scale + scalescale^2/2),1:2000,standat$L,standat$U,standat$D)
  yhate <- rowMeans(yhate)
  
  plusinc_configs <- unique(standat$X[indi,])
  etai <- matrix(alpha,nrow=2000,ncol=nrow(plusinc_configs)) + beta %*% t(plusinc_configs)
  yhati <- array(dim = dim(etai))
  for (k in 1:nrow(plusinc_configs))
    yhati[,k] <- Yhatify(etai[,k],exp(scale + scalescale^2/2),1:2000,standat$L,standat$U,standat$D)
  yhati <- rowMeans(yhati)
  
  marginaldiff <- rbind(marginaldiff,cbind(evbal=i,post_summary_dt(yhate-yhat0),type="Exculpatory",source="Linear model"))
  marginaldiff <- rbind(marginaldiff,cbind(evbal=i,post_summary_dt(yhati-yhat0),type="Inculpatory",source="Linear model"))
  
}
marginaldiff <- marginaldiff[!((evbal==-3 & type == "Exculpatory") | (evbal==3 & type == "Inculpatory"))]


#### Experiment 2 ####
whichtask <- "exculpatory_conditional"
source("~/code/casereveal/analysis/process_data.R")
ldat_cond <- makelegaldat(scendat,subjdat,clickdat)
setkey(ldat_cond,uid,question)

rfitint_cond <- list()
rfit_cond <- list()
for (c in evcond) {
  standat <- makestandat(ldat_cond[cond_evidence==c],interact = T)
  rfit_cond[[c]] <- sampling(model,standat,chains=4,init=transfunc_init,iter=1000,warmup=500,
                                pars=c(basepars,rpars))
  rfitint_cond[[c]] <- sampling(intmodel,standat,chains=4,init=transfunc_init,iter=1000,warmup=500,
                                   pars=c(basepars,rpars,intpars))
}

bfit_cond <- list()
for (c in evcond) {
  standat <- makestandat(ldat_cond[cond_evidence==c],binresp = T)
  bfit_cond[[c]] <- sampling(bardmodel,standat,chains=4,iter=1000,warmup=500,pars=basepars)
}

##### Time trend ####
bdat_cond <- merge(ldat_cond,makebalancedat(ldat_cond),by=c("uid","scenario"))
fit_trend <- lmer(rating ~ 1 + cond_evidence*I(n_inculp-n_exculp)*question +
                    (1 + I(n_inculp-n_exculp)*question | uid) + (1|scenario) ,data=bdat_cond,REML=F)

##### Scenario-level differences ####
expose_stan_functions(rfit_cond$balanced)
beta_scen_resp_diff <- get_beta_scen_resp(rfit_cond$balanced) - get_beta_scen_resp(rfit_cond$credible)
alpha_scen_resp_diff <- get_alpha_scen_resp(rfit_cond$credible) - get_alpha_scen_resp(rfit_cond$balanced)

##### learning model ####
rldat <- makestanrldat(ldat_cond)
rlmodel <- stan_model("/home/seth/code/legalmodels/refcase_transfunc1_flat.stan")
rlfit <- optimizing(rlmodel,rldat,as_vector=F)

#### Experiment 3 ####

whichtask <- "exculpatory_rateless"
source("~/code/casereveal/analysis/process_data.R")
ldat <- makelegaldat(scendat,subjdat,clickdat)

bfit_rate <- list(list(),list())
names(bfit_rate) <- ratecond
for (ec in evcond) {
  for (rc in ratecond) {
    standat <- makestandat(ldat[cond_evidence==ec & cond_rating==rc],binresp = T)
    bfit_rate[[rc]][[ec]] <- sampling(bardmodel,standat,chains=4,iter=1000,warmup=500,pars=basepars)
  }
}

rfit_rate <- list()
for (c in evcond) {
  standat <- makestandat(ldat[cond_evidence==c & cond_rating=="with"])
  rfit_rate[[c]] <- sampling(model,standat,chains=4,init=transfunc_init,iter=1000,warmup=500,
                             pars=c(basepars,rpars))
}
