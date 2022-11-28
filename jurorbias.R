whichtask <- "exculpatory_burdenofproof"
source("~/code/casereveal/analysis/prep_stan.R")
# qualdat <- fread(paste0(fpath, "_qualtrics.csv"))
# qualdat <- qualdat[uid %in% dat$uid]
# jbscale <- qualdat[,str_subset(names(qualdat),"AC\\d+")[1:12], with=F]
# jbscale <- sapply(jbscale,as.numeric)
# 
# revcode <- c(11,12)
# jbscale[,revcode] <- jbscale[,revcode]*-1 + 6
# for (i in 1:ncol(jbscale)) {
#   ji <- jbscale[,i]
#   ji[is.na(ji)] <- mean(ji, na.rm = T)
#   jbscale[,i] <- ji
# }
# jbdat <- merge(dat, data.table(uid = qualdat$uid, bias = scale(rowMeans(jbscale))), by="uid")

#visualize
#####
evidat <- dat[,.(n_exculp=str_detect(sapply(.(physical,document,witness,character),as.character),"ex") %>% sum(),
                 n_inculp=str_detect(sapply(.(physical,document,witness,character),as.character),"in") %>% sum(),
                 rating,bias,legalguilt),by=.(uid,scenario)]
evidat[, c("evidence_balance","biasgroup") := 
         .(n_inculp - n_exculp, cut_number(bias, 3, ordered_result=T, labels=c("low","medium","high")))]
ggplot(evidat,aes(x=evidence_balance,y=rating, color=biasgroup)) + 
  geom_smooth(method="lm", se=F, size=0.5) + 
  geom_point(data=evidat[,.(rating=mean(rating),N=.N),by=.(evidence_balance,biasgroup)],
             aes(size=N),position = position_dodge(.25), shape=20) +
  scale_color_discrete("juror bias") + scale_size(breaks=c(10,150))

evidat[,.(evidence_balance = seq(-4,4,.1), guilty = glm(legalguilt ~ evidence_balance, family = binomial) %>%
            predict(newdata = list(evidence_balance = seq(-4,4,.1)), type = "response")),
       by=biasgroup] %>% ggplot(aes(y=guilty,x=evidence_balance,color=biasgroup)) + geom_line() +
  geom_point(data=evidat[,.(guilty=mean(legalguilt),N=.N),by=.(evidence_balance,biasgroup)],
             aes(size=N),position = position_dodge(.25), shape=20) + scale_size(breaks=c(10,150)) +
  scale_color_discrete("juror bias") + scale_y_continuous("% voting guilty",labels = scales::percent)

#guilt judgements
#####
library(lme4)
dat[, severity := scenario/31]
fit <- glmer(guilty ~ 1 + standard + standard*severity + standard*bias + 
               standard*physical + standard*document + standard*witness + standard*character +
               (1 + standard | uid) + (1 + standard | scenario), 
      family = binomial, data = guiltdat)

bfit <- glmer(legalguilt ~ 1 + bias + 
                bias*physical + bias*document + bias*witness + bias*character +
                (1 | uid) + (1 | scenario), 
              family = binomial, data = dat[uid %in% r2pass])

#####
Z <- model.matrix(~ 1 + physical + document + witness + character + bias +
                    bias:physical + bias:document + bias:witness + bias:character, data = dat)
Z <- Z[,-c(1:ncol(standat_base$X))]
standat <- c(standat_base, D=1, P0=ncol(Z), Z = list(Z))
model <- stan_model("~/code/legalmodels/models/mixtransfer_uncsigma_fixedeff.stan")
rtfit <- sampling(model,standat,chains=4,init=0,iter=1000,warmup=500,thin=2,
                pars=c("alpha","beta_mu","sigma_subj","sigma_scen","beta_subj","beta_scen","sigma","eta"))

model <- stan_model("~/code/legalmodels/models/logistic_fixedeff.stan")
standat <- c(standat_base, P0=ncol(Z), Z = list(Z), Y = list(dat[,as.numeric(subjguilt)]))
ljfit <- sampling(model,standat,chains=4,init=0,iter=1000,warmup=500,thin=2,
                 pars=c("alpha","beta_mu","sigma_subj","sigma_scen","beta_subj","beta_scen","eta"))

model <- stan_model("~/code/legalmodels/models/logistic_mixedfactor.stan")
standat <- c(standat_base, Y=list(dat$legalguilt), K=1, B=list(dat[,bias[1],by=uid][,V1]))
foo <- sampling(model,standat,chains=4,iter=2000,warmup=1000,thin=2,
                pars=c("Lambda","sigma_lambda","beta_mu",
                       "sigma_subj_resid","sigma_scen","log_lik"))

standat$K <- 2
foo2 <- sampling(model,standat,chains=4,iter=2000,warmup=1000,thin=2,
                pars=c("Lambda","sigma_lambda","beta_mu",
                       "sigma_subj_resid","sigma_scen","log_lik"))

##### 
juh <- stan_plot(foo,"Lambda")$data
# juh <- stan_plot(foo2,str_subset(names(foo2),"Lambda\\[\\d+,2\\]"))$data
juh$type <- c("juror_bias","baseline",rep(c("physical","document","witness","character"),c(3,3,3,2)))
juh$type <- ordered(juh$type,levels=unique(juh$t))
juh$level <- c("NA","NA",c(rep(c("clear_ex","ambiguous","clear_in"),3),c("clear_ex","clear_in")))
juh$level <- ordered(juh$level,c(evidord,"NA"))
ggplot(juh,aes(y=mean,x=type,fill=level)) + 
  geom_pointrange(aes(ymin=ll,ymax=hh),position = position_dodge(width=0.3), shape=21) +
  xlab(NULL) + ylab("factor loading") + geom_hline(yintercept = 0) + 
  theme(legend.title = element_blank())
stan_hist(foo,"r2")
#####
bias_ev_plt <- function(fit) {
  samps <- extract(fit,c("beta_mu","alpha"))
  highbias <- lowbias <- array(dim=dim(samps$beta_mu))
  for (i in 1:1000) {
    highbias[i,] <- samps$beta_mu[i,] %*% diag(12) + samps$alpha[i,] %*% diag(12)
    lowbias[i,] <- samps$beta_mu[i,] %*% diag(12) + samps$alpha[i,] %*% (diag(12)*-1)
  }
  colnames(highbias) <- colnames(standat$X)
  colnames(lowbias) <- colnames(standat$X)
  highbias <- melt(as.data.table(highbias))[,.(effect=mean(value),
                                               lb=quantile(value,0.025),
                                               ub=quantile(value,0.975)),by=variable]
  lowbias <- melt(as.data.table(lowbias))[,.(effect=mean(value),
                                             lb=quantile(value,0.025),
                                             ub=quantile(value,0.975)),by=variable]
  biasdat <- rbind(cbind(highbias,juror_bias="high"),cbind(lowbias,juror_bias="low"))
  biasdat$type <- rep(c("baseline","physical","document","witness","character"),c(1,3,3,3,2)) %>% rep(2)
  biasdat$level <- c("NA",rep(c("clear_ex","ambiguous","clear_in"),3),c("clear_ex","clear_in")) %>% rep(2)
  biasdat[,c("level","juror_bias"):=.(ordered(level,c(evidord,"NA")), ordered(juror_bias,levels=c("low","medium","high")))]
  return(biasdat)
}

biasdat <- rbind(bias_ev_plt(ljfit)[,measure:="rating"],bias_ev_plt(ljfit)[,measure:="vote"])
ggplot(biasdat[type!="baseline"],aes(y=effect,x=type,fill=juror_bias)) + 
  geom_pointrange(aes(ymin=lb,ymax=ub),position = position_dodge(width=0.3), shape=21) + 
  xlab(NULL) + geom_hline(yintercept = 0) + facet_wrap(c("measure","level"), scales = "free_y") + scale_fill_discrete("juror bias")
