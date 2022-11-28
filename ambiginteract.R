whichtask <- "exculpatory_burdenofproof"
source("~/code/casereveal/analysis/process_data.R")
source("~/code/casereveal/analysis/prep_stan.R")
#####
#X <- X[,!str_detect(colnames(X), "character")]
# standat <- c(standat_base, Pe=3, Pi=3, Pa=3,
#              Y=list(dat[,as.numeric(legalguilt)]),
#              Xe=list(X[,colnames(X[,str_detect(colnames(X), "_ex")])]),
#              Xi=list(X[,colnames(X[,str_detect(colnames(X), "_in")])]),
#              Xa=list(X[,colnames(X[,str_detect(colnames(X), "ambiguous")])]))
# model <- stan_model("~/code/legalmodels/hailmary.stan")
# wtf <- sampling(model,standat,chains=4,iter=2500,warmup=1000,thin=2,
#                 pars=c("mu_w_exc","mu_w_inc",
#                        "mu_beta_0","mu_beta_inc","mu_beta_exc","mu_beta_amb",
#                        "sigma_0","sigma_inc","sigma_exc","sigma_amb",
#                        "mu_beta_axe","mu_beta_axi","sigma_axe","sigma_axi",
#                        "mu_eps","sigma_eps"))

newsubj <- subjdat[uid %in% subjdat[start>as.POSIXct("2020-03-01",format="%Y-%m-%d"),uid],unique(uid)]
subjsets <- list(subjdat[!(uid %in% newsubj),unique(uid)],newsubj,subjdat[,unique(uid)])

#####
#full interaction model
fit_allint <- list()
for (i in 1:length(subjsets)) {
  standat <- makestandat(makelegaldat(scendat,subjdat,subset = subjsets[[i]]),interact = T, binresp = T)
  model <- stan_model("~/code/legalmodels/logistic_interactions_lapse.stan")
  pars <- c("beta_mu","sigma_subj","sigma_scen","sigma_lambda_mu","sigma_lambda_subj",
            "lambda_mu","mu_eps","sigma_eps")
  if (i==3) pars <- c(pars,"log_lik")
  fit_allint[[i]] <- sampling(model,standat,chains=4,iter=2000,warmup=500,thin=2,pars=pars)
}
model <- stan_model("~/code/legalmodels/logistic_lapse.stan")
standat <- makestandat(makelegaldat(scendat,subjdat),interact = F, binresp = T)
fit_null <- sampling(model,standat,chains=4,iter=2000,warmup=500,thin=2,
                     pars=c("beta_mu","sigma_subj","sigma_scen","mu_eps","sigma_eps","log_lik"))
#####
#visualize in raw data
evplt <- data.table()
for (i in 1:length(subjsets)) {
  evidat <- dat[uid %in% subjsets[[i]],.(n_exculp=str_detect(sapply(.(physical,document,witness,character),as.character),"ex") %>% sum(),
                                         n_inculp=str_detect(sapply(.(physical,document,witness,character),as.character),"in") %>% sum(),
                                         n_ambig=str_detect(sapply(.(physical,document,witness,character),as.character),"amb") %>% sum(),
                                         phys_amb=(physical=="ambiguous"),doc_amb=(document=="ambiguous"),wit_amb=(witness=="ambiguous"),
                                         rating,legalguilt),by=.(uid,scenario)]
  evidat[,evbal:=n_inculp-n_exculp]
  evplt <- rbind(evplt,
                 evidat[,.(evbal = seq(-4,4,.1), guilty = glm(legalguilt ~ evbal, family = binomial) %>%
                             predict(newdata = list(evbal = seq(-4,4,.1)), type = "response"),
                           subj=c("old subjects","new subjects","all subjects")[i]),
                        by=n_ambig>0])
}
ggplot(evplt,aes(y=guilty,x=evbal,color=n_ambig)) + geom_line() +
  geom_point(data=evidat[,.(guilty=mean(legalguilt),N=.N),by=.(evbal,n_ambig>0)],
             aes(size=N),position = position_dodge(.25), shape=20) + facet_wrap("subj",nrow=1)

#####
#minimal interaction model
X <- model.matrix(~ 1 + n_exculp + n_inculp + phys_amb + doc_amb + wit_amb +
                    n_exculp:phys_amb + n_exculp:doc_amb + n_exculp:wit_amb + 
                    n_inculp:phys_amb + n_inculp:doc_amb + n_inculp:wit_amb, data=evidat)
P <- dim(X)[2]
standat <- list(Nsubj=Nsubj, Nscen=Nscen, N=N, P=P, Scen=C, Subj=S, X=X, Y=evidat[,as.numeric(legalguilt)])
model <- stan_model("~/code/legalmodels/logistic_lapse.stan")
fit_disag_amb <- sampling(model,standat,chains=4,iter=2000,warmup=500,
                pars=c("beta_mu","sigma_subj","sigma_scen","mu_eps","sigma_eps"))
effdat <- as.data.table(stan_plot(fit_disag_amb,"beta_mu")$data)
tmp <- str_split_fixed(colnames(X),":",2)
effdat[,c("direction","type") := .(tmp[,1], tmp[,2])]

X <- model.matrix(~ 1 + n_exculp + n_inculp + n_ambig + n_ambig:n_exculp + n_ambig:n_inculp, data=evidat)
P <- dim(X)[2]
standat <- list(Nsubj=Nsubj, Nscen=Nscen, N=N, P=P, Scen=C, Subj=S, X=X, Y=evidat[,as.numeric(legalguilt)])
model <- stan_model("~/code/legalmodels/logistic_lapse.stan")
fit_agg_amb <- sampling(model,standat,chains=4,iter=2000,warmup=500,
                          pars=c("beta_mu","sigma_subj","sigma_scen","mu_eps","sigma_eps"))

effdat2 <- as.data.table(stan_plot(fit_agg_amb,"beta_mu")$data)
tmp <- str_split_fixed(colnames(X),":",2)
effdat2[,c("direction","type") := .(tmp[,1], tmp[,2])]
effdat3 <- rbind(effdat,effdat2)[type!=""][,type:=ordered(type,levels=unique(type))]

ggplot(effdat3,aes(y=mean,x=type,color=direction)) + 
  geom_pointrange(aes(ymin=l,ymax=h),position = position_dodge(0.5),size=0.75) + 
  geom_errorbar(aes(ymin=ll,ymax=hh),width=0,size=0.25,position = position_dodge(0.5)) + 
  geom_hline(yintercept =0) + ylab("interaction effect size") + scale_color_discrete("interaction with...") +
  scale_x_discrete("ambiguous evidence type",labels=c("physical","document","witness","total # ambiguous"))
