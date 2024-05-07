library(tidybayes)
whichtask <- "exculpatory_conditional"
source("miscfunctions.R")
source("prep_stan.R")
source("process_data.R")
basepars <- c("mu_alpha","mu_beta","sigma_alpha_subj","sigma_alpha_scen","sigma_beta_subj","sigma_beta_scen",
              "alpha_scen","beta_scen","alpha_subj","beta_subj","log_lik")
rpars <- c("mu_alpha_resp","mu_beta_resp","mu_scale","sigma_scale","scale_subj","Yhat")
condmodel <- stan_model("models/transfunc1_conditions.stan")
ldat_cond <- makelegaldat(scendat,subjdat,clickdat)
setkey(ldat_cond,uid,question)

ldat_cond[,evconf:=paste0(physical,document,witness,character)]
ldat_cond_cred <- ldat_cond[evconf %in% ldat_cond[cond_evidence=="credible",evconf]]

standat <- makestandat(ldat_cond_cred,cond="cond_evidence",levelref=T)
rfit_cond_cred <- sampling(condmodel,standat,chains=4,iter=500,warmup=200,pars=c(basepars,rpars))
rfit_cond_cred <- recover_types(rfit_cond_cred,standat$ref)

ldat_cond[,evconf:=paste0(physical,document,witness,character)]
ldat_cond_def <- ldat_cond[evconf %in% ldat_cond[cond_evidence=="defenseless",evconf]]

standat <- makestandat(ldat_cond_def,cond="cond_evidence",levelref=T)
rfit_cond_def <- sampling(condmodel,standat,chains=4,iter=500,warmup=200,pars=c(basepars,rpars))
rfit_cond_def <- recover_types(rfit_cond_def,standat$ref)

foo <- rbind(
  cbind(maineff_cond(rfit_cond),fit_on="All cases"),
  cbind(maineff_cond(rfit_cond_cred),fit_on="Credible only"),
  cbind(maineff_cond(rfit_cond_def)[valence %in% c("Baseline","Inculpatory")],fit_on="Defenseless only")
)

false_god <- ggplot(foo[fit_on=="All"],aes(y=mean,ymin=lb,ymax=ub,x=cond,color=cond)) + 
  geom_pointrange(position = position_dodge(width=0.25)) + facet_wrap(vars(valence), scales="free", nrow=1) + 
  scale_color_manual("Condition",values=condscheme) + scale_fill_manual("Condition",values=condscheme) + scale_x_discrete(NULL, labels=NULL) +
  ylab("Case strength (points)")

harbinger <- ggplot(foo[cond=="Balanced"],aes(y=mean,ymin=lb,ymax=ub,x=fit_on,color=cond)) + geom_pointrange() + 
  facet_wrap(vars(valence), scales="free", nrow=1) + scale_color_manual("Condition",values=condscheme) + 
  scale_fill_manual("Condition",values=condscheme) + ylab("Case strength (points)") +  xlab(NULL) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.8, hjust = 0.8))

cast_down <- ggplot(foo,aes(y=mean,ymin=lb,ymax=ub,x=fit_on,color=cond)) + geom_pointrange(position = position_dodge(width=0.5)) + 
  facet_wrap(vars(valence), scales="free", nrow=1) + scale_color_manual("Condition",values=condscheme) + scale_fill_manual("Condition",values=condscheme) + 
  ylab("Case strength (points)") +  xlab(NULL) + theme(axis.text.x = element_text(angle = 45, vjust = 0.8, hjust = 0.8)) 
