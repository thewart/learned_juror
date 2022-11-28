# whichtask <- "exculpatory_burdenofproof"
source("~/code/casereveal/analysis/prep_stan.R")
ns <- 1e3
initf <- function() return(list(sigma=runif(1,15,30)))
yppars1 <- c("sigma","theta")
yppars2 <- c("w_trans","l_trans","s_trans")

model <- stan_model("~/code/legalmodels/models/mixtransfer_uncsigma.stan")
standat <- c(standat_base, D=1)
fit_1sigma_1k_cdflik <- sampling(model,standat,chains=4,init=0,iter=1500,warmup=500,thin=2,
                pars=c("mu","eta","tau","sigma","theta","log_lik"))

model <- stan_model("~/code/legalmodels/models/sv_n.stan")
fit_linear <- sampling(model,standat_base,chains=4,init=initf,iter=1500,warmup=500,thin=2,
                       pars=c("mu","eta","tau","sigma","theta","log_lik"))

model <- stan_model("~/code/legalmodels/models/mixtransfer_2sigma.stan")
expose_stan_functions(model)
standat <- c(standat_base, K=1, ns=ns)
fit_2sigma_1k <- sampling(model,standat,chains=4,init=initf,iter=3500,warmup=500,thin=3,
                          pars=c("mu","eta","tau","sigma","w_trans","l_trans","s_trans","theta","log_lik"))
standat$K <- 2
fit_2sigma_2k <- sampling(model,standat,chains=4,init=initf,iter=3500,warmup=500,thin=3,
                          pars=c("mu","eta","tau","sigma","w_trans","l_trans","s_trans","theta","log_lik"))

model <- stan_model("~/code/legalmodels/models/mixtransfer_1sigma.stan")
standat <- c(standat_base, K=1, ns=ns)
fit_1sigma_1k <- sampling(model,standat,chains=4,init=initf,iter=1500,warmup=500,thin=2,
                          pars=c("mu","eta","tau","sigma","w_trans","l_trans","s_trans","theta","log_lik"))
standat$K <- 2
fit_1sigma_2k <- sampling(model,standat,chains=1,init=initf,iter=5000,warmup=1000,thin=2,
                          pars=c("mu","eta","tau","sigma","w_trans","l_trans","s_trans","theta","log_lik"),
                          control=list(max_treedepth=12))

baseline <- waic(extract_log_lik(fit_linear))
wdat <- rbind(compare(baseline,waic(extract_log_lik(fit_1sigma_1k))),
              compare(baseline,waic(extract_log_lik(fit_1sigma_2k))),
              compare(baseline,waic(extract_log_lik(fit_2sigma_1k))),
              compare(baseline,waic(extract_log_lik(fit_2sigma_2k)))) %>% as.data.table()
wdat$noise <- c("external","external","external & internal","external & internal")
wdat$mixture <- c("one component","two components","one component","two components")
ggplot(wdat,aes(y=elpd_diff,x=noise,color=mixture)) + 
  geom_pointrange(aes(ymin=elpd_diff-2*se,ymax=elpd_diff+2*se),position=position_dodge(width=0.25))

ppdat <- rbind(data.table(rating=dat$rating, model="empirical"),
               data.table(rating=as.vector(rypred(extract(fit_linear,pars=yppars1),"linear")),model="linear"),
               data.table(rating=as.vector(rypred(extract(fit_2sigma_1k,pars=c(yppars1,yppars2)),"doublenoise")),model="internal & external noise"),
               data.table(rating=as.vector(rypred(extract(fit_1sigma_2k,pars=c(yppars1,yppars2)),"externalnoise")),model="external noise"),
               data.table(rating=as.vector(rypred(extract(fit_1sigma_1k_cdflik,pars=yppars1),"internalnoise")),model="internal noise"))
ggplot(ppdat[model %in% c("empirical","linear","internal noise")],aes(x=rating,y=stat(density))) + geom_histogram(bins=101) + facet_wrap(facets = "model", ncol = 1)

qplot(summary(fit_linear,"theta")[[1]][,1])



model <- stan_model("~/code/legalmodels/models/betacdf_lik.stan")
fit_beta <- sampling(model,standat,chains=1,init=0,iter=200,warmup=100,thin=2,
                                 pars=c("mu","eta","tau","inv_scale","theta","log_lik"))
