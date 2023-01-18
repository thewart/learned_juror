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
