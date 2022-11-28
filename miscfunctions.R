printci <- function(dat,lev,d=2,n=F,wrap=F) {
  r <- dat[level==lev,1:3,with=F]
  if (n) {
    r <- lapply(r,function(x) -1*x)
    CI <- paste0("CI=[", round(r$ub,d), ", ", round(r$lb,d),"]")
  } else {
    CI <- paste0("CI=[", round(r$lb,d), ", ", round(r$ub,d),"]")
  }
  
  if (wrap) {
    out <- paste0(round(r$mean,d)," (",CI,")")
  } else {
    out <- paste0("mean=",round(r$mean,d),", ",CI)
  }
  return(out)
}

printci2 <- function(r,d=2,n=F,wrap=F) {
  if (n) {
    r <- -r
    CI <- paste0("CI=[", round(ub(r),d), ", ", round(lb(r),d),"]")
  } else {
    CI <- paste0("CI=[", round(lb(r),d), ", ", round(ub(r),d),"]")
  }
  
  if (wrap) {
    out <- paste0(round(mean(r),d)," (",CI,")")
  } else {
    out <- paste0("mean=",round(mean(r),d),", ",CI)
  }
  return(out)
}


poisbinom <- function(lambda, theta, n = 1e6) {
  N <- rpois(n,lambda)
  out <- rnk(N,theta,n)
  out$lambda <- lambda
  return(out)
}

unifbinom <- function(a, b, theta, n = 1e6) {
  N <- rcategorical(n,rep(1/(b-a),b-a+1)) + a - 1
  out <- rnk(N,theta,n)
  out$a <- a
  out$b <- b
  return(out)
}

normbinom <- function(mu, sigma, theta, n = 1e6) {
  N <- rnorm(n,mu,sigma) %>% round()
  N[N<0] <- 0
  out <- rnk(N,theta,n)
  out$mu <- mu
  out$sigma <- sigma
  return(out)
}

rcategorical <- function(n,prob) {
  return(rmultinom(n,1,prob) %>% apply(2,function(x) which(x==1)))
}

idunno <- function(n,k,p) {
  factorial(n)/factorial(n-k) * p^n
}

rnk <- function(N, theta, n=1e6) {
  K <- rbinom(n,N,theta)
  out <- data.table(N=N,K=K,theta=theta)
  setkey(out,K)
  return(out)
}

mvwnft <- function(data,ws,cond=NA) {
  mvplt <- data.table()
  for (i in 0:(29-ws)) {
    mvfit <- lmer(rating ~ I(n_inculp - n_exculp) +  (1 + I(n_inculp - n_exculp) | uid),
                  data=data[question %in% c(i:(i+ws-1))])
    mvfix <- summary(mvfit)$coefficients
    iplt <- data.table(question=median(c(i:(i+ws-1))),
                                    bl=mvfix[1,1],bl_se=mvfix[1,2],
                                    evbal=mvfix[2,1],evbal_se=mvfix[2,2])
    # if (nrow(mvfix)==3) {
    #   iplt <- cbind(iplt,exceff=mvfix[3,1],exceff_se=mvfix[3,2])
    # } else {
    #   iplt <- cbind(iplt,exceff=NA,exceff_se=NA)
    # }
    mvplt <- rbind(mvplt,iplt)
  }
  mvplt$condition <- cond
  return(mvplt)
}


post_summary_dt <- function(samps,name = colnames(samps)) {
  if (length(dim(samps)) %in% c(0,1))
    samps <- matrix(samps,ncol=1)
  return(data.table(mean=colMeans(samps),lb=apply(samps,2,quantile,prob=0.025),
                    ub=apply(samps,2,quantile,prob=0.975),level=name))
}


aggsumm_cond <- function(fit,ordcond,condlabel="condition",label=NULL,cont=NULL,resp_scale=F) {
  dt <- data.table()
  for (cond in ordcond) {
    dti <- aggsumm(fit[[cond]],resp_scale=resp_scale) %>% post_summary_dt()
    dti[[condlabel]] <- cond
    dt <- rbind(dt,dti)
  }
  if (!is.null(label)) dt[[label]] <- cont
  return(dt)
}

aggsumm <- function(fit,incind=c(3,6,9,11),excind=c(1,4,7,10),ambind=c(2,5,8),resp_scale=F) {
  if (resp_scale) {
    beta <- extract(fit,pars="mu_beta_resp")[[1]]
    alpha <- extract(fit,pars="mu_alpha_resp")[[1]]
  } else {
    beta <- extract(fit,pars="mu_beta")[[1]]
    alpha <- extract(fit,pars="mu_alpha")[[1]]
  }
  
  samps <- data.table(baseline=alpha,exculpatory=rowMeans(beta[,excind]),ambiguous=rowMeans(beta[,ambind]),
                 inculpatory=rowMeans(beta[,incind]),combined=rowMeans(beta))
  return(samps)
}

aggsummlabels <- function(dt, condlevels=c("Balanced","Credible","Defenseless"),
                          evlevels=c("Baseline","Exculpatory","Ambiguous","Inculpatory","Combined")) {
  dt$condition <- factor(tools::toTitleCase(dt$condition), levels=condlevels)
  dt$level <- factor(tools::toTitleCase(dt$level),levels = evlevels)
  return(dt)
}

makestanrldat <- function(ldat, rescale=F) {
  standat <- makestandat(ldat)
  standat$Trial <- ldat$question + 1
  standat$Cond <- ldat[,ordered(cond_evidence) %>% as.numeric()]
  standat$NX <- rowSums(standat$X)
  if (rescale) standat$Y <- scale(standat$Y) %>% as.vector()
  return(standat)
}

combine_ab <- function(fit,p1="mu_alpha",p2="mu_beta") {
  return(cbind(extract(fit,p1)[[1]],extract(fit,p2)[[1]]))
}

label_dt <- function(dt) {
  evtypes <- c("Baseline","Physical","Document","Witness","Character")
  evlevels <- c("Baseline","Exculpatory","Ambiguous","Inculpatory")
  dt$type <- rep(evtypes,times=c(1,3,3,3,2)) %>% factor(levels=evtypes)
  dt$level <- c(evlevels[1],rep(evlevels[-1],4)[-11]) %>% factor(levels=evlevels)
  return(dt)
}

transfunc_init <- function(priorscale=0.25,P=11) {
  return(list(mu_alpha=rnorm(1),mu_beta=rnorm(P),
              sigma_alpha_scen=abs(rnorm(1,0,priorscale)),sigma_alpha_subj=abs(rnorm(1,0,priorscale)),
              sigma_beta_scen=abs(rnorm(P,0,priorscale)),sigma_beta_subj=abs(rnorm(P,0,priorscale)),
              sigma_mu_lambda=abs(rnorm(1,0,0.05)),sigma_lambda_subj=abs(rnorm(1,0,0.05)),
              mu_scale=rnorm(1,0,0.5),sigma_scale=abs(rnorm(1,0,.25))))
}

# postpredsamp <- function(X, Subj, Scen, mu_alpha, mu_beta, alpha_subj, alpha_scen, beta_subj, beta_scen) {
#   beta_subj <- beta_subj %>% purrr::array_branch(1)
#   beta_scen <- beta_scen %>% purrr::array_branch(1)
#   
#   return(etaize(X, Subj, Scen, mu_alpha, alpha_subj, alpha_scen, mu_beta, beta_subj, beta_scen))
# }

resample <- function(bl, subj=T, scen=F, scale=T, binint=F) {
  nsamp <- length(bl$mu_alpha)
  if (subj) {
    bl$alpha_subj <- array(dim=c(nsamp,bl$Nsubj))
    bl$beta_subj <- array(dim=c(nsamp,bl$Nsubj,bl$P))
    if (binint) bl$lambda_subj <- array(dim=c(nsamp,bl$Nsubj,bl$P2))
  }
  if (scen) {
    bl$alpha_scen <- array(dim=c(nsamp,bl$Nscen))
    bl$beta_scen <- array(dim=c(nsamp,bl$Nscen,bl$P))
  }
  if (scale) bl$scale_subj <- array(dim=c(nsamp,bl$Nsubj))
  
  for (i in 1:nsamp) {
    if (subj) {
      bl$alpha_subj[i,] <-  rnorm(bl$Nsubj)*bl$sigma_alpha_subj[i]
      bl$beta_subj[i,,] <- diag(bl$sigma_beta_subj[i,]) %*% (rnorm(bl$Nsubj*bl$P) %>% matrix(nrow=bl$P))
      if (binint) bl$lambda_subj[i,,] <- bl$sigma_lambda_subj[i]*rnorm(bl$Nsubj*bl$P2) %>% matrix(nrow=bl$P2)
    }
    if (scen) {
      bl$alpha_scen[i,] <-  rnorm(bl$Nscen)*bl$sigma_alpha_scen[i]
      bl$beta_scen[i,,] <- diag(bl$sigma_beta_scen[i,]) %*% (rnorm(bl$Nscen*bl$P) %>% matrix(nrow=bl$P))
    }
    if (scale) bl$scale_subj[i,] <- exp(bl$mu_scale[i] + bl$sigma_scale[i] * rnorm(bl$Nsubj))
  }
  
  return(postpredeta(bl,binint,T))
}

postpredeta <- function(bl, binint=F, list_out=F) {
  nsamp <- length(bl$mu_alpha)
  bl$eta <- matrix(nrow=nsamp,ncol=bl$N)
  
  for (i in 1:nsamp) {
    arglist <- with(bl,list(X, Subj, Scen, mu_alpha[i], alpha_subj[i,], alpha_scen[i,], mu_beta[i,],
                            beta_subj[i,,] %>% purrr::array_branch(1), beta_scen[i,,] %>% purrr::array_branch(1)))
    if (binint) arglist <- c(arglist[1],list(bl$Z),arglist[-1],list(bl$mu_lambda[i,],bl$lambda_subj[i,,] %>% purrr::array_branch(1)))
    bl$eta[i,] <- do.call(etaize,arglist)
  }
  
  if (list_out) return(bl)
  else return(bl$eta)
}

postpredsample <- function(bl, rfunc=rlogis, list_out=F, resample_subj=F, resample_scen=F, resample_scale=resample_subj, binint=F) {
  nsamp <- length(bl$mu_alpha)
  if (!("eta" %in% names(bl)))
    bl <- resample(bl, resample_subj, resample_scen, resample_scale, binint) %>% postpredeta(binint, T)
  if (!("eps" %in% names(bl))) 
    bl$eps <- rfunc(nsamp*bl$N) %>% matrix(nrow=nsamp,ncol=bl$N)
  
  bl$yrep <- matrix(nrow=nsamp,ncol=bl$N)
  for (i in 1:nsamp)  bl$yrep[i,] <- pnorm(bl$eta[i,] + bl$eps[i,], 0, bl$scale_subj[i,bl$Subj])
  bl$yrep <- (bl$yrep*(bl$U - bl$L) + bl$L) %>% round()
  
  if (list_out) return(bl)
  else return(bl$yrep)
}

postpredyhat <- function(bl, list_out=F) {
  nsamp <- length(bl$mu_alpha)
  bl$yhat <- matrix(nrow=nsamp,ncol=bl$N)
  for (i in 1:nsamp) bl$yhat[i,] <- with(bl,Yhatify(eta[i,],scale_subj[i,],Subj,L,U,D))
  
  if (list_out) return(bl)
  else return(bl$yhat)
}

postpredeta_cfact <- function(bl, comp, list_out=F, binint=F, rfunc=rlogis) {
  nsamp <- length(bl$mu_alpha)
  for (c in comp) bl[[c]] <- bl[[c]] * 0
  bl$eta <- postpredeta(bl,binint)
  # bl <- postpredsample(bl, rfunc=rfunc, binint=binint, list_out = T)
  
  if (list_out) return(bl)
  else return(bl$eta)
}

ppc_ydist <- function(y,yrep,binsize) {
  nbins=100/binsize
  yrepbin <- cut(yrep,seq(0,100,binsize),1:nbins,include.lowest = T) %>% as.numeric() %>% array(dim=dim(yrep))
  yrepcount <- sapply(1:nbins,function(x) rowSums(yrepbin==x))
  yrepdt <- data.table(rating=1:nbins,proportion=colMeans(yrepcount),
                       lb=apply(yrepcount,2,quantile,prob=0.025),ub=apply(yrepcount,2,quantile,prob=0.975))
  ydt <- data.table(rating=1:nbins,
                    proportion=sapply(1:nbins,function(x) sum((cut(y,seq(0,100,binsize),1:nbins,include.lowest = T) %>% as.numeric())==x)))
  return(list(y=ydt,yrep=yrepdt))
}

get_beta_scen_resp <- function(fit) {
  scalescale <- extract(fit,"sigma_scale")[[1]]
  scale <- extract(fit,"mu_scale")[[1]]
  beta_scen <- extract(fit,"beta_scen")[[1]]
  beta <- extract(fit,"mu_beta")[[1]]
  for (i in 1:nrow(beta)) beta_scen[i,,] <- sweep(beta_scen[i,,],2,beta[i,],FUN = "+")
  for (i in 1:dim(beta_scen)[2])
    for (j in 1:dim(beta_scen)[3])
      beta_scen[,i,j] <- Yhatify(beta_scen[,i,j],exp(scale + scalescale^2/2),1:2000,0,100,1)
  return(beta_scen)
}

get_alpha_scen_resp <- function(fit) {
  scalescale <- extract(fit,"sigma_scale")[[1]]
  scale <- extract(fit,"mu_scale")[[1]]
  alpha_scen <- sweep(extract(fit,"alpha_scen")[[1]],1,extract(fit,"mu_alpha")[[1]],FUN = "+")
  for (i in 1:ncol(alpha_scen)) alpha_scen[,i] <- Yhatify(alpha_scen[,i],exp(scale + scalescale^2/2),1:2000,0,100,1)
  return(alpha_scen)
}

intvsmain <- function(standat,fitint) {
  intconts <- str_split_fixed(colnames(standat$Z),":",n=2) %>% apply(1,function(y) colnames(standat$X) %in% y)
  mainsum <- extract(fitint,"mu_beta")[[1]] %*% intconts %>% post_summary_dt()
  colnames(foo) <- c("compsum","clb","cub")
  return(cbind(foo,extract(fitint,"mu_lambda")[[1]] %>% post_summary_dt))
}

calc_con_ind <- function(I,B) {
  X <- rep(0,length(B))
  X[I] <- 1
  return(calc_con(X,B))
}
calc_se_ind <- function(I,S) {
  X <- rep(0,nrow(S))
  X[I] <- 1
  return(calc_con_se(X,S))
}
calc_t_ind <- function(I,B,S) {
  X <- rep(0,length(B))
  X[I] <- 1
  return(calc_con(X,B)/calc_con_se(X,S))
}

calc_con <- function(X,B) (t(B) %*% X) |> as.vector()
calc_con_se <- function(X,S) (t(X) %*% S %*% X) |> diag() |> sqrt()
