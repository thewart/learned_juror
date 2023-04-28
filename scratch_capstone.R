whichtask <- "capstone_prolific"
source("process_data.R")
source("miscfunctions.R")
source("prep_stan.R")
ldat_cap <- makelegaldat(scendat,subjdat,clickdat)
model <- stan_model("models/transfunc1_intervention_conditions.stan")
# model2 <- stan_model("models/transfunc1_intervention_interaction.stan")

standat <- makestandat(ldat_cap,cond="cond_evidence",levelref=T)
standat$Z <- ldat_cap[,cond_capstone]

rfit_cap <- sampling(model,standat,iter=1000,warmup=500,chains=4)
rfit_cap <- recover_types(rfit_cap,standat$ref)

standat <- makestandat(ldat[cond_evidence=="balanced"])
standat$Z <- ldat[cond_evidence=="balanced",cond_capstone]
fit_bal <- sampling(model,standat,iter=500,warmup=200,chains=4,init=transfunc_init_intervention)
# standat$W <- data.table(standat$X %*% summary(rfit,"mu_beta")[[1]][,1],standat$Subj)[,sum(V1),by=V2]$V1 |> scale() |> as.vector()

standat <- makestandat(ldat[cond_evidence=="credible"])
standat$Z <- ldat[cond_evidence=="credible",cond_capstone]
fit_cred <- sampling(model,standat,iter=500,warmup=200,chains=4,init=transfunc_init_intervention)

samps_bal <- extract(rfit_cap$balanced,c("mu_lambda","mu_alpha","mu_beta","mu_alpha_resp","mu_beta_resp","mu_scale","sigma_scale"))
samps_cred <- extract(rfit_cap$credible,c("mu_lambda","mu_alpha","mu_beta","mu_alpha_resp","mu_beta_resp","mu_scale","sigma_scale"))


instruct_dat <- rbind(data.table(Baseline=samps_bal$mu_alpha_resp,
                                 Evidence=rowMeans(samps_bal$mu_beta_resp)) |> 
                        post_summary_dt() |> cbind(instructed=F,context="Balanced"),
                      data.table(Baseline=with(samps_bal,pointscale(mu_alpha+mu_lambda[,1],mu_scale,sigma_scale)),
                                 Evidence=with(samps_bal,apply(sweep(mu_beta,1,-mu_lambda[,2]),2,pointscale,mu_scale,sigma_scale,T)) |> rowMeans()) |> 
                        post_summary_dt() |> cbind(instructed=T,context="Balanced"),
                      data.table(Baseline=samps_cred$mu_alpha_resp,
                                 Evidence=rowMeans(samps_cred$mu_beta_resp)) |> 
                        post_summary_dt() |> cbind(instructed=F,context="Credible"),
                      data.table(Baseline=with(samps_cred,pointscale(mu_alpha+mu_lambda[,1],mu_scale,sigma_scale)),
                                 Evidence=with(samps_cred,apply(sweep(mu_beta,1,-mu_lambda[,2]),2,pointscale,mu_scale,sigma_scale,T)) |> rowMeans()) |> 
                        post_summary_dt() |> cbind(instructed=T,context="Credible"))

ggplot(instruct_dat,aes(y=mean,ymin=lb,ymax=ub,x=instructed,color=context,group=context)) + geom_pointrange(position = position_dodge(width=.5)) + 
  facet_wrap(vars(level),scales="free",strip.position = "bottom") + geom_line(position = position_dodge(width=.5)) + 
  scale_x_discrete(NULL,labels=c("Pre-instruction","Post-instruction")) + ylab("Case strength (points)") + theme(strip.placement = "outside", legend.position = "top")

qcheck <- function(whichtask) {
  source("process_data.R",local = T)
  ldat <- makelegaldat(scendat,subjdat,clickdat)
  bdat <- merge(ldat,makebalancedat(ldat),by=c("uid","scenario"))
  if (whichtask=="exculpatory_rateless") bdat <- bdat[cond_rating=="with"]
  out <- bdat[,.(tstat=cor.test(rating,I(n_inculp - n_exculp))$statistic),by=uid]
  out$task <- whichtask
  return(out)
}

foo <- rbind(qcheck("capstone"),qcheck("exculpatory_burdenofproof"),qcheck("exculpatory_conditional"),qcheck("exculpatory_rateless"))
ggplot(foo,aes(x=tstat,color=task)) + geom_density()


quickndirty <- function(whichtask,evc=c("credible","random"),cpc=0:1,thresh=-Inf)  {
  source("process_data.R",local = T)
  ldat <- makelegaldat(scendat,subjdat,clickdat)
  bdat <- merge(ldat,makebalancedat(ldat),by=c("uid","scenario"))
  tdat <- bdat[,.(tstat=cor.test(rating,I(n_inculp - n_exculp))$statistic),by=uid]
  fit <- lmer(rating ~ 1 + n_inculp + n_exculp + n_ambig + (1 + n_inculp + n_exculp + n_ambig||uid),
       data=bdat[cond_evidence %in% evc & cond_capstone %in% cpc & uid %in% tdat[tstat>thresh,uid]]) |> summary()
  return(cbind(as.data.table(fit$coefficients), valence=c("baseline","inculpatory","exculpatory","ambiguous"),
               task=whichtask,thresh=thresh,cond_evidence=evc,cond_capstone=cpc))
}

threshes <- c(-Inf,0,1,2)
evcond <- c("credible","balanced")
foo <- data.table()
for (t in threshes)
  for (c in evcond) 
    foo <- rbind(foo,quickndirty(whichtask,c,c(0,1),t))
foo[,valence:=ordered(valence,levels=unique(valence))]
foo[,notbaseline:=valence!="baseline"]
ggplot(foo, aes(x=valence,y=Estimate,color=cond_evidence)) +
  geom_pointrange(aes(ymin=Estimate - `Std. Error`,ymax=Estimate + `Std. Error`),position = position_dodge(width=0.3)) + 
  geom_hline(data=data.table(y=c(0,NA),notbaseline=c(T,F)),aes(yintercept=y)) + xlab(NULL) + scale_color_brewer(NULL,palette = "Dark2") + scale_y_continuous("Case strength (points)",breaks=seq(-10,80,20)) + 
  facet_grid(notbaseline ~ thresh,scales="free_y",drop=F,labeller = "label_both") + theme(axis.text.x=element_text(angle=30,hjust=1),strip.text.y = element_blank(),panel.spacing = unit(1,"lines"))

juh <- data.table()
for (t in c(1,1.5,2))
  for (c in c(0,1))
    juh <- rbind(juh,quickndirty(whichtask,evcond,c,t))
juh[,valence:=ordered(valence,levels=unique(valence))]
juh[,notbaseline:=valence!="baseline"]
ggplot(juh, aes(x=valence,y=Estimate,color=as.logical(cond_capstone))) +
  geom_pointrange(aes(ymin=Estimate - `Std. Error`,ymax=Estimate + `Std. Error`),position = position_dodge(width=0.3)) + geom_hline(data=data.table(y=c(0,NA),notbaseline=c(T,F)),aes(yintercept=y)) +
  xlab(NULL) + scale_color_brewer(NULL,palette = "Dark2") + scale_y_continuous("Case strength (points)",breaks=seq(-10,80,20)) + 
  facet_grid(notbaseline ~ thresh,scales="free_y",drop=F,labeller = "label_both") + theme(axis.text.x=element_text(angle=30,hjust=1),strip.text.y = element_blank(),panel.spacing = unit(1,"lines"))
