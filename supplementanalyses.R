##### PPC -- Case strength rating distribution ####

expose_stan_functions(rfit)
ob <- c(standat,extract(rfit))
yrep <- postpredsample(ob)

yrepfreq <- sapply(0:100,function(x) rowMeans(yrep==x))
yrepdt <- data.table(rating=0:100,proportion=colMeans(yrepfreq),
                    lb=apply(yrepfreq,2,quantile,prob=0.025),ub=apply(yrepfreq,2,quantile,prob=0.975))
ydt <- data.table(rating=0:100,proportion=sapply(0:100,function(x) mean(ldat$rating==x)))
ggplot(yrepdt,aes(x=rating,y=proportion)) + geom_ribbon(alpha=0.33,aes(ymin=lb,ymax=ub)) + geom_line(linetype=2) +
  geom_line(data=ydt) + ylab("Frequency") + xlab("Rating")

binsize <- 5
ind <- c(72,107,7)
ppcdt <- list(y=data.table(),yrep=data.table())
# ind <- sample.int(standat$Nsubj,6)
for (i in ind) {
  iind <- standat$Subj==i
  ippcdt <- ppc_ydist(ldat$rating[iind],yrep[,iind],binsize)
  ippcdt$y$Subj <- paste("Subject",i)
  ippcdt$yrep$Subj <- paste("Subject",i)
  ppcdt$y <- rbind(ppcdt$y,ippcdt$y)
  ppcdt$yrep <- rbind(ppcdt$yrep,ippcdt$yrep)
}
ggplot(ppcdt$yrep,aes(x=rating,y=proportion)) + geom_line(linetype=2) +
  geom_col(data=ppcdt$y,alpha=0.5) + ylab("Frequency") + xlab("Rating") + facet_wrap("Subj")


yrepbin <- cut(yrep, breaks = seq(0,100,10), include.lowest = T, labels = F) |> matrix(nrow = nrow(yrep))
calibdat <- data.table()
for (i in 1:nrow(yrep)) calibdat <- data.table(Y=ldat$rating,Yrep=yrepbin[i,])[,.(Yhat=median(Y)),by=Yrep] |> cbind(iter=i) |> rbind(calibdat)
calibdat[,post_summary_dt(Yhat),by=Yrep] |> ggplot(aes(y=mean,ymin=lb,ymax=ub,x=(Yrep-.5)*10)) + geom_pointrange() + geom_abline()

##### interactions
vardecomp <- (1-apply(cyhat,1,var)/apply(extract(rfitint,"Yhat")[[1]][],1,var))
cor(colMeans(cyhat),summary(rfit,"Yhat")[[1]][,1])
qplot(y=colMeans(cyhat),x=summary(rfitint,"Yhat")[[1]][,1],alpha=I(0.1)) + geom_abline() +
  xlab("Rating (full model)") + ylab("Rating (without interactions)")

intvsmainplt <- function(standat,fitint) {
  intconts <- str_split_fixed(colnames(standat$Z),":",n=2) %>% apply(1,function(y) colnames(standat$X) %in% y)
  foo <- extract(fitint,"mu_beta")[[1]] %*% intconts %>% post_summary_dt()
  colnames(foo) <- c("compsum","clb","cub")
  return(cbind(foo,extract(fitint,"mu_lambda")[[1]] %>% post_summary_dt))
}
intdt <- intvsmainplt(standat,rfitint)
ggplot(intdt,aes(x=compsum,y=mean)) + geom_point() +
  geom_errorbar(aes(ymin=lb,ymax=ub),alpha=0.25,width=0) + geom_errorbarh(aes(xmin=clb,xmax=cub),alpha=0.25,height=0) +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + ylab("Interaction effect (a.u.)") + xlab("Sum of main effects (a.u.)")

##### really seth
standat <- makestandat(ldat)
powmodel_fixed <- stan_model("models/transfunc1_pow_fixed.stan")
pfit <- sampling(powmodel_fixed,standat,chains=4,iter=500,warmup=250,pars=c("mu_alpha","mu_beta","mu_scale","sigma_scale","mu_gamma","Yhat","log_lik"),
                 control=list(max_treedepth=14))


##### maaath ####
math_xlab <- "Observed evidence (*M*)"
math_ylab <- "Inferred evidence (*N-M*|*M*)"
tm <- 4
q <- .5
ff <- c(1/3,2/3,3/2,3)
cmp_dat <- data.table()
for (f in ff) {
  juh <- rcmpbinom2(tm,tm*f,q,1e6)[,.(khat=mean(n-m),.N),by=m]
  juh$ff <- f
  cmp_dat <- rbind(cmp_dat,juh)
}
cmp_dat <- cmp_dat[m < cmp_dat[N<500,min(m)],-"N",with=F]
cmp_dat <- rbind(cmp_dat,data.table(m=0:max(cmp_dat$m),khat=tm*(1-q),ff=1))
cmp_dat$dist <- "Conway-Maxwell Poisson"
cmp_plt <- ggplot(cmp_dat,aes(y=khat,x=m,color=as.ordered(ff))) + geom_point() + geom_line() + coord_cartesian(c(0,NA),c(0,NA)) + ylab(math_ylab) + xlab(math_xlab) +
  scale_color_ordinal(NULL,guide=guide_legend(ncol = 2),labels=c("Strong under","Weak under","Poisson disperson","Weak over-disperson","Strong over")) + ggtitle("CMP evidence collection") +
  theme(legend.position = "top", legend.justification = "right")

plotpq <- function(p,q,nu=6) {
  theta <- (p*(1-q))/(1-p*q)
  return(data.table(khat=(nu-0:nu)*theta,m=0:nu,theta=paste0("p=",p,",q=",q)))
}
binom_dat <- rbind(plotpq(1,0),plotpq(.5,.5),plotpq(0,1)) |> cbind(dist="Binomial")
binom_plt <-  ggplot(binom_dat,aes(y=khat,x=m,color=as.ordered(theta) |> forcats::fct_rev())) + geom_line() + facet_wrap(vars(dist)) +
  scale_color_discrete(NULL,labels=c("p=1, q<1","p=q=0.5", "p<1, q=1")) + scale_x_continuous(math_xlab, labels=c("0",NA,NA,expression(nu))) + 
  scale_y_continuous(math_ylab, labels=c("0",NA,NA,expression(nu))) + theme(axis.title.x = ggtext::element_markdown(), axis.title.y = ggtext::element_markdown())


nbkhat <- function(mu,ff,q) {
  v = mu*ff
  nu = mu^2/(v-mu)
  theta = (1-q)*mu/(nu+q*mu)
  return(data.table(m=0:6,khat=(nu+0:6)*theta,ff=ff))
}
nbdat <- rbind(nbkhat(4,3,.5),
      nbkhat(4,3/2,.5),
      data.table(m=0:6,khat=rep(2,7),ff=1)) |> 
  cbind(dist="Negative binomial")

unbound_dat <- rbind(nbdat,cmp_dat)
unbound_dat[,dist:=factor(dist,levels=unique(dist))]
unbound_plt <- ggplot(unbound_dat[m<7],aes(y=khat,x=m,color=as.ordered(ff))) + geom_point() + geom_line() + facet_wrap(vars(dist),nrow=1) + 
  scale_color_ordinal("Dispersion",labels=c("Strong under","Weak under","Normal","Weak over","Strong over")) +
  coord_cartesian(c(0,NA),c(0,NA)) + ylab(math_ylab) + xlab(math_xlab) + theme(axis.title.x = ggtext::element_markdown(), axis.title.y = ggtext::element_markdown())

plot_grid(binom_plt,unbound_plt,ncol=2,rel_widths = c(0.6,1),axis="tb",align = "h",labels = "auto")

