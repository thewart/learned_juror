whichtask <- "exculpatory_burdenofproof"
source("process_data.R")
ldat <- makelegaldat(scendat,subjdat,clickdat)
setkey(ldat,uid,question)
bdat <- merge(ldat,makebalancedat(ldat),by=c("uid","scenario"))

foo <- stan_lmer(rating ~ 1 + n_inculp + n_exculp + n_ambig + 
                   (1 + n_inculp + n_exculp + n_ambig | uid) + 
                   (1 + n_inculp + n_exculp + n_ambig || scenario), 
                 data = bdat, iter=400, warmup=200, chains=4)

fixed <- as.matrix(foo, pars=c("(Intercept)","n_inculp","n_exculp","n_ambig"))
baseline <- fixed[,1]
evav <- rowMeans(fixed[,2:4])

scen_dat <- data.table()
for (i in 1:31) {
  scen_samps <- as.matrix(foo,regex_pars=paste0("\\:",i,"\\]"))
  scen_evav <- evav + rowMeans(scen_samps[,2:4])
  scen_baseline <- baseline + scen_samps[,1]
  scen_dat <- rbind(scen_dat, post_summary_dt(cbind(scen_evav,scen_baseline),ci=.9) |> cbind(scenario=i))
}

scen_dat$dname <- fread("../task/data/Scenarios.csv", select = "Defendant name")[[1]] |> rep(each=2)
scen_dat[,dname := ordered(dname,levels=unique(dname))]
scen_dat$ethnic <- F
scen_dat[dname %in% c("Diego Alvarez", "Juan Rodriguez", "Sergei Ivankov", "Mohammed Farooqi"), ethnic := T]
ggplot(scen_dat,aes(y=mean,ymin=lb,ymax=ub,x=dname,color=ethnic)) + geom_pointrange() + facet_wrap(vars(level),scales = "free_y",ncol=1) + 
  scale_x_discrete(NULL, guide = guide_axis(angle=50)) + theme_cowplot() + ylab("Effect size (points)") + theme(legend.position = "none")

subj_dat <- data.table()
for (i in ldat[,unique(uid)]) {
  subj_samps <- as.matrix(foo,regex_pars=i)
  subj_evav <- evav + rowMeans(subj_samps[,2:4])
  subj_baseline <- baseline + subj_samps[,1]
  subj_dat <- rbind(subj_dat, post_summary_dt(cbind(subj_evav,subj_baseline),ci=.9) |> cbind(uid=i))
}
