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
