printci <- function(dat,lev,d=2,n=F,wrap=F) {
  r <- dat[valence==lev,.(mean,lb,ub)]
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
  out <- rnm(N,theta,n)
  out$lambda <- lambda
  return(out)
}

unifbinom <- function(a, b, theta, n = 1e6) {
  N <- rcategorical(n,rep(1/(b-a),b-a+1)) + a - 1
  out <- rnm(N,theta,n)
  out$a <- a
  out$b <- b
  return(out)
}

normbinom <- function(mu, sigma, theta, n = 1e6) {
  N <- rnorm(n,mu,sigma) %>% round()
  N[N<0] <- 0
  out <- rnm(N,theta,n)
  out$mu <- mu
  out$sigma <- sigma
  return(out)
}

rnbinombinom <- function(nu,p,q,nsim = 1e6) {
  n <- rnbinom(nsim,nu,1-p)
  out <- rnm(n,q,nsim)
  out$nu <- nu
  out$p <- p
  return(out)
}

rcmpbinom <- function(lambda,nu,q,nsim = 1e6) {
  n <- rcmp(nsim,lambda,nu)
  out <- rnm(n,q,nsim)
  out$nu <- nu
  out$lambda <- lambda
  return(out)
}

CMPare <- function(theta,target_mean,target_var) 
  (target_mean - ecmp(theta[1],theta[2]))^2 + (target_var - vcmp(theta[1],theta[2]))^2

rcmpbinom2 <- function(mean,var,q,nsim = 1e6) {
  theta <- optim(c(5,1),CMPare,target_mean=mean,target_var=var)$par
  return(rcmpbinom(theta[1],theta[2],q,nsim))
}

dnbinombinom_analytic <- function(x,nu,m,p,q) {
  pm <- (1-q)*p
  return(dnbinom(x,nu+m,1-pm))
}

enbinombinom <- function(nu,m,p,q) {
  pm <- (1-q)*p
  return((nu+m)*pm/(1-pm) + m)
}

rlogbinom <- function(p,q,nsim=1e6) {
  n <- rlogarithmic(nsim,p)
  out <- rnm(n,q,nsim)
  out$p <- p
  return(out)
}

rcategorical <- function(n,prob) {
  return(rmultinom(n,1,prob) %>% apply(2,function(x) which(x==1)))
}

idunno <- function(n,k,p) {
  factorial(n)/factorial(n-k) * p^n
}

rnm <- function(n, q, nsim=1e6) {
  m <- rbinom(nsim,n,q)
  out <- data.table(n=n,m=m,q=q)
  setkey(out,m)
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

post_summary_dt <- function(samps,name = colnames(samps), ci=.95) {
  if (length(dim(samps)) %in% c(0,1))
    samps <- matrix(samps,ncol=1)
  return(data.table(mean=colMeans(samps),lb=apply(samps,2,quantile,prob=(1-ci)/2),
                    ub=apply(samps,2,quantile,prob=0.5+ci/2),level=name))
}

aggsumm_cond <- function(fit,ordcond,condlabel="condition",label=NULL,cont=NULL,resp_scale=F) {
  dt <- data.table()
  for (cond in ordcond) {
    dti <- aggsumm1(fit[[cond]],resp_scale=resp_scale) %>% post_summary_dt()
    dti[[condlabel]] <- cond
    dt <- rbind(dt,dti)
  }
  if (!is.null(label)) dt[[label]] <- cont
  return(dt)
}

aggsumm1 <- function(fit,incind=c(3,6,9,11),excind=c(1,4,7,10),ambind=c(2,5,8),resp_scale=F) {
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

aggsumm2 <- function(fit,incind=c(3,6,9,11),excind=c(1,4,7,10),ambind=c(2,5,8),resp_scale=F) {
  if (resp_scale) {
    #   beta <- extract(fit,pars="mu_beta_resp")[[1]]
    alpha <- extract(fit,pars="mu_alpha_resp")[[1]]
  } else {
    alpha <- extract(fit,pars="mu_alpha")[[1]]
  }
  beta <- extract(fit,pars="mu_beta")[[1]]
  
  samps <- data.table(exculpatory=rowMeans(beta[,excind]),ambiguous=rowMeans(beta[,ambind]),inculpatory=rowMeans(beta[,incind]),combined=rowMeans(beta))
  
  if (resp_scale) {
    scale <- extract(fit,pars=c("mu_scale","sigma_scale"))
    samps <- lapply(samps,pointscale,scale$mu_scale,scale$sigma_scale,center=T) |> as.data.table()
  }
  
  samps <- cbind(baseline=alpha,samps)
  return(samps)
}

aggsummlabels <- function(dt, condlevels=c("Balanced","Credible","Defenseless"),
                          evlevels=c("Baseline","Exculpatory","Ambiguous","Inculpatory","Combined")) {
  dt$condition <- factor(tools::toTitleCase(dt$condition), levels=condlevels)
  dt$level <- factor(tools::toTitleCase(dt$level),levels = evlevels)
  return(dt)
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
  nsamp <- length(scale)
  for (i in 1:nrow(beta)) beta_scen[i,,] <- sweep(beta_scen[i,,],2,beta[i,],FUN = "+")
  for (i in 1:dim(beta_scen)[2])
    for (j in 1:dim(beta_scen)[3])
      beta_scen[,i,j] <- Yhatify(beta_scen[,i,j],exp(scale + scalescale^2/2),1:nsamp,0,100,1)
  return(beta_scen)
}

get_alpha_scen_resp <- function(fit) {
  scalescale <- extract(fit,"sigma_scale")[[1]]
  scale <- extract(fit,"mu_scale")[[1]]
  nsamp <- length(scale)
  alpha_scen <- sweep(extract(fit,"alpha_scen")[[1]],1,extract(fit,"mu_alpha")[[1]],FUN = "+")
  for (i in 1:ncol(alpha_scen)) alpha_scen[,i] <- Yhatify(alpha_scen[,i],exp(scale + scalescale^2/2),1:nsamp,0,100,1)
  return(alpha_scen)
}

get_alpha_subj_resp <- function(fit) {
  scale <- extract(fit,"scale_subj")[[1]]
  alpha_subj <- sweep(extract(fit,"alpha_subj")[[1]],1,extract(fit,"mu_alpha")[[1]],FUN = "+")
  nsamp <- nrow(alpha_subj)
  for (i in 1:ncol(alpha_subj)) alpha_subj[,i] <- Yhatify(alpha_subj[,i],scale[,i],1:nsamp,0,100,1)
  return(alpha_subj)
}

get_beta_subj_resp <- function(fit) {
  samps <- extract(fit,c("mu_beta","beta_subj"))
  scale <- extract(fit,"scale_subj")[[1]]
  beta_subj <- sweep(samps$beta_subj,c(1,3),samps$mu_beta,"+")
  nsamp <- nrow(beta_subj)
  for (i in 1:dim(beta_subj)[2])
    for (j in 1:dim(beta_subj)[3])
      beta_subj[,i,j] <- Yhatify(beta_subj[,i,j],scale[,i],1:nsamp,0,100,1)
  return(beta_subj-50)
}

get_beta_subj_bias <- function(fit) {
  samps <- extract(fit,c("mu_beta","beta_subj"))
  beta_subj <- sweep(samps$beta_subj,c(1,3),samps$mu_beta,"+")
  return(apply(beta_subj,c(1,2),sum)/apply(abs(beta_subj),c(1,2),sum))
}

get_beta_subj_tot <- function(fit) {
  samps <- extract(fit,c("mu_beta","beta_subj"))
  beta_subj <- sweep(samps$beta_subj,c(1,3),samps$mu_beta,"+")
  return(apply(abs(beta_subj),c(1,2),sum))
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

pointscale <- function(x,mu_scale,sigma_scale,center=F) {
  y <- pnorm(x/exp(mu_scale+sigma_scale^2/2))*100
  if (center) y <- y - 50
  return(y)
}

factor_inorder <- function(x) return(factor(x,levels=unique(x)))

parse_evidence <- function(X,baseline=T) {
  X[str_detect(evidence,"physical"),type:="Physical"]
  X[str_detect(evidence,"document"),type:="Document"]
  X[str_detect(evidence,"witness"),type:="Witness"]
  X[str_detect(evidence,"character"),type:="Character"]
  X[str_detect(evidence,"clear_ex"),valence:="Exculpatory"]
  X[str_detect(evidence,"ambiguous"),valence:="Ambiguous"]
  X[str_detect(evidence,"clear_in"),valence:="Inculpatory"]
  # X[,c(".chain",".iteration","evidence",".variable"):=NULL]
  if (baseline) X[is.na(type),c("type","valence"):="Baseline"]
  
  X[,type:=factor(type,levels=c("Baseline","Physical","Document","Witness","Character"))]
  X[,valence:=factor(valence,levels=c("Baseline","Exculpatory","Ambiguous","Inculpatory"))]
  return(X)
}

parse_capstone <- function(X) {
  X[str_detect(.variable,"alpha"),type:="Baseline"]
  X[str_detect(.variable,"beta"),type:="Ave. evidence"]
  X[str_detect(.variable,"pre"),capped:=F]
  X[str_detect(.variable,"post"),capped:=T]
  
  return(X)
}

weights_by_cond_plot <- function(X,colorscale,bard=F) {
  ytitle <- "Case strength (points)"
  if (bard) ytitle <- "Guilty judgment (log odds)"
  plt <- ggplot(X, aes(x=valence,y=mean,color=cond)) + geom_pointrange(aes(ymin=lb,ymax=ub),position = position_dodge(width=0.3)) + 
    geom_hline(data=data.table(y=c(0,NA),notbaseline=c(T,F)),aes(yintercept=y)) +
    xlab(NULL) + scale_color_manual(NULL,breaks=names(colorscale),values=colorscale) + scale_y_continuous(ytitle) + theme(axis.text.x=element_text(angle=30,hjust=1))
}

weights_within_valence <- function(X) return(X[,mean(.value),by=.(cond,valence,.draw)][,post_summary_dt(V1),by=.(cond,valence)])
average_weights <- function(X,summary=T) {
  out <- X[valence!="Baseline",.(.value=mean(.value),valence="Ave. evidence"),by=.(cond,.draw)]
  if (summary) {
    return(out[,post_summary_dt(.value),by=.(cond,valence)])
  } else {
    return(out)
  }
}

set_factor_context <- function(x) return(factor(str_to_title(x),levels = c("Balanced","Credible","Defenseless")))
cull_defenseless <- function(dt) return(dt[!(cond=="Defenseless" & valence %in% c("Exculpatory","Ambiguous","Ave. evidence"))])
set_factor_rate <- function(DT) {
  x <- DT$cond
  context_cond <- str_split_i(x,":",1) |> set_factor_context()
  DT[,cond:=context_cond]
  context_rate <- str_split_i(x,":",2)
  DT[,rating:=factor(str_to_title(context_rate),levels=c("Without","With"))]
  return(DT)
}

set_factor_valence <- function(x) return(factor(str_to_title(x),levels = c("Baseline","Exculpatory","Ambiguous","Inculpatory")))

subjeff_cond <- function(stanfit,cond_df,resp_scale=T,average_over="type",post_int=F) {
  subj_samps <- gather_draws(stanfit,alpha_subj[uid],beta_subj[uid,evidence]) |> setDT() |> parse_evidence()
  mu_samps <- gather_draws(stanfit,mu_alpha[cond],mu_beta[cond,evidence]) |> setDT() |> parse_evidence()
  
  subj_samps <- subj_samps[cond_df,on="uid"]
  subj_samps <- subj_samps[mu_samps,on=c(".draw","cond","type","valence")]
  subj_samps[,.value:= .value + i..value]
  
  if (average_over=="type") {
    subj_samps <- subj_samps[,mean(.value),by=.(uid,cond,valence,.draw)]
  } else if (average_over=="all") {
    subj_samps <- subj_samps[,mean(.value),by=.(uid,cond,valence=="Baseline",.draw)]
    subj_samps[,valence:=factor(valence,levels=c(TRUE,FALSE),labels=c("Baseline","Ave. evidence"))]
  }
  
  if (post_int) {
    int_samps <- gather_draws(stanfit,lambda_subj[uid,i]) |> setDT()
    muint_samps <- gather_draws(stanfit,mu_lambda[cond,i]) |> setDT()
    
    int_samps <- int_samps[cond_df,on="uid"]
    int_samps <- int_samps[muint_samps,on=c(".draw","cond","i")]
    int_samps[,.value:= .value + i..value]
    
    subj_samps <- rbind(
      subj_samps[valence=="Baseline"][int_samps[i==1,.(uid,.draw,int=.value)], on=c("uid",".draw")],
      subj_samps[valence!="Baseline"][int_samps[i==2,.(uid,.draw,int=.value)], on=c("uid",".draw")]
    )
    subj_samps[,V1 := V1 + int]
    subj_samps[,int := NULL]
  }
  
  if (resp_scale) {
    scale_samps <- gather_draws(stanfit,scale_subj[uid]) |> setDT()
    subj_samps <- subj_samps[scale_samps[,.(uid,.draw,scale=.value)], on=c("uid",".draw")]
    subj_samps[,V1:=pnorm(V1/scale)*100 - 50*(valence!="Baseline")]
  }

  subj_eff <- subj_samps[,.(cond=cond[1],mean=mean(V1)),by=.(uid,valence)]
  
  subj_eff[,cond:=set_factor_context(cond)]
  return(cull_defenseless(subj_eff))
}

maineff_cond <- function(stanfit,resp_scale=T,average_over="type",post_int=F) {
  mu_samps <- gather_draws(stanfit,mu_alpha[cond],mu_beta[cond,evidence]) |> setDT() |> parse_evidence()
  
  if (average_over=="type") {
    mu_samps <- mu_samps[,mean(.value),by=.(cond,valence,.draw)]
  } else if (average_over=="all") {
    mu_samps <- mu_samps[,mean(.value),by=.(cond,valence=="Baseline",.draw)]
    mu_samps[,valence:=factor(valence,levels=c(TRUE,FALSE),labels=c("Baseline","Ave. evidence"))]
  }
  
  if (post_int) {
    int_samps <- gather_draws(stanfit,mu_lambda[cond,i]) |> setDT()
    
    mu_samps <- rbind(
      mu_samps[valence=="Baseline"][int_samps[i==1,.(cond,.draw,int=.value)], on=c("cond",".draw")],
      mu_samps[valence!="Baseline"][int_samps[i==2,.(cond,.draw,int=.value)], on=c("cond",".draw")]
    )
    mu_samps[,V1 := V1 + int]
    mu_samps[,int := NULL]
  }
  
  if (resp_scale) {
    scale_samps <- spread_draws(stanfit,mu_scale,sigma_scale) |> setDT()
    mu_samps <- mu_samps[scale_samps[,.(.draw,mu_scale,sigma_scale)],on=".draw"]
    mu_samps[,.value:=pnorm(V1/exp(mu_scale + sigma_scale^2/2))*100 - 50*(valence!="Baseline")]
  }
  
  mu_samps[,cond:=set_factor_context(cond)]
  return(mu_samps[,post_summary_dt(.value),by=.(cond,valence)] |> cull_defenseless())
}
