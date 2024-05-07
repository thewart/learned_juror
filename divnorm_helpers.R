ybarify <- function(bdat, resp="binary") {
  ybar_dat <- bdat[,.(ybar=ifelse(resp=="binary", mean(bardguilt), mean(rating)/100),
                      se=ifelse(resp=="binary", sd(bardguilt)/sqrt(.N), sd(rating/100)/sqrt(.N)),
                      M=.N), by=.(physical,document,witness,character,n_exculp,n_inculp,n_ambig)]
  ybar_dat[,evconf:=paste0(physical,document,witness,character)]
  ybar_dat[,nev:= n_inculp + n_exculp + n_ambig]
  return(ybar_dat)
}

fit_onesided <- function(ybar_dat, model) {
  ybar_onesided <- get_onesided(ybar_dat)
  X <- model.matrix(~ physical + document + witness + character, ybar_onesided)[,-1]
  X_pred <- model.matrix(~ physical + document + witness + character, ybar_dat)[,-1]
  standat <- list(N=nrow(X), P=ncol(X), X=X, ybar=ybar_onesided$ybar, M=ybar_onesided$M, m0=1, se=ybar_onesided$se,
                  X_pred=X_pred, N_pred=nrow(X_pred))
  return(optimizing(model, standat, as_vector=F, hessian=T))
}

fit_singleton <- function(ybar_dat, model) {
  ybar_single <- ybar_dat[nev < 2]
  X <- model.matrix(~ physical + document + witness + character, ybar_single)[,-1]
  X_pred <- model.matrix(~ physical + document + witness + character, ybar_dat)[,-1]
  standat <- list(N=nrow(X), P=ncol(X), X=X, ybar=ybar_single$ybar, M=ybar_single$M, se=ybar_single$se, m0=1,
                  X_pred=X_pred, N_pred=nrow(X_pred))
  return(optimizing(model, standat, as_vector=F))
}

get_onesided <- function(df) return(df[((n_exculp>0) + (n_inculp>0) + (n_ambig>0)) %in% c(0,1)])

plot_ypred_poly <- function(ybar_dat, ypred, degree=3) {
  ybar_comp <- cbind(ybar_dat, yhat=ypred)
  plt <- ggplot(ybar_comp[nev>0], aes(y=ybar, x=yhat, color=ordered(nev))) + geom_point(alpha=0.3) + geom_abline() + 
    geom_smooth(method="lm", formula=y ~ poly(x, degree), se=F) + xlab("Sum of evidence weights (compressed)") + ylab("Observed case strength") +
    scale_colour_viridis_d("Evidence \n count", end=0.9)
  return(plt)
}

plot_ypred_poly_onesided <- function(ybar_dat, ypred, degree=3, evidscheme=NULL) {
  if (is.null(evidscheme)) evidscheme <- c("Baseline"="black","Inculpatory"="#ba4040ff","Exculpatory"="#406bbaff", "Ambiguous"="#765884ff")
  ybar_dat <- valence_from_evcounts(ybar_dat)
  ybar_comp <- cbind(ybar_dat, yhat=ypred)
  plt <- ggplot(ybar_comp, aes(y=ybar, x=yhat)) + geom_pointrange(aes(ymin=ybar-2*se, ymax=ybar+2*se, color=valence)) + geom_abline(alpha=0.5, linetype=2) +
    geom_smooth(method="lm", formula= y ~ poly(x, 3), se=F, color="black") + scale_color_manual("Evidence valence", breaks=names(evidscheme)[-1], values=evidscheme) +
    ylab("Observed case strength") + xlab("Sum of weights")
  return(plt)
}

lm_est_int <- function(ybar_dat, ypred, degree=1) {
  ybar_comp <- cbind(ybar_dat, yhat=ypred)
  return(ybar_comp[nev>0,lm(ybar-yhat ~ 1 + poly(yhat, degree) + poly(yhat, degree):nev)])
}

lm_est_resid <- function(ybar_dat, ypred, degree=1) {
  ybar_comp <- cbind(ybar_dat, yhat=ypred)
  return(ybar_comp[nev>0,lm(ybar-yhat ~ 1 + poly(yhat, degree))])
}

test_int <- function(ybar_dat, ypred, degree=1) {
  ybar_comp <- cbind(ybar_dat, yhat=ypred)
  return(anova(lm_est_resid(ybar_dat, ypred, degree), lm_est_int(ybar_dat, ypred, degree)))
}

fit_full <- function(ybar_dat, model) {
  X <- X_pred <- model.matrix(~ physical + document + witness + character, ybar_dat)[,-1]
  standat <- list(N=nrow(X), P=ncol(X), X=X, ybar=ybar_dat$ybar, M=ybar_dat$M, m0=1, se=ybar_dat$se,
                  X_pred=X_pred, N_pred=nrow(X_pred))
  return(optimizing(model, standat, as_vector=F, hessian=T))
}

fit_divglm <- function(ybar_dat, model, type="count") {
  X <- model.matrix(~ physical + document + witness + character, ybar_dat)[,-1]
  if (type=="count") {
    Z <- model.matrix(~ factor(nev), data=ybar_dat)[, -1]
  } else if (type=="flex") {
    Z <- X
  } else if (type=="linear_count") {
    Z <- model.matrix(~ nev, data=ybar_dat)[, -1]
    dim(Z) <- c(length(Z), 1)
  }
  standat <- list(N=nrow(X), P=ncol(X), X=X, ybar=ybar_dat$ybar, M=ybar_dat$M, m0=1, se=ybar_dat$se, Z=Z, Q=ncol(Z))
  return(optimizing(model, standat, as_vector=F, hessian=T))
}

valence_from_evcounts <- function(dt) {
  dt[n_exculp>0, valence := "Exculpatory"]
  dt[n_inculp>0, valence := "Inculpatory"]
  dt[n_ambig>0, valence := "Ambiguous"]
  dt[is.na(valence), valence := "Baseline"]
  dt[, valence := factor(valence, levels=c("Baseline","Exculpatory","Ambiguous","Inculpatory"))]
  return(dt)
}

compare_models <- function(model1, model2, parind1=4, parind2=3) {
  kdiff <- purrr::list_c(model1$par[1:parind1]) |> length() - purrr:::list_c(model2$par[1:parind2]) |> length()
  lldiff <- model1$value - model2$value
  lrtest <- 1-pchisq(2*lldiff, kdiff)
  lr_pred <- exp(lldiff - kdiff) 
  return(data.table(lldiff=lldiff, kdiff=kdiff, pval=lrtest, lr_pred=lr_pred))
}

pdiv <- function(q, pars, nev) return(pt(q/c(exp(pars$gamma*nev)), df=exp(pars$lognu)))

maketable <- function(ybar_dat, shrinky_model, divnorm_model, response) {
  shrinky_fit <- fit_full(ybar_dat, shrinky_model)
  div_count_fit <- fit_divglm(ybar_dat, divnorm_model, type="linear_count")
  div_flex_fit <- fit_divglm(ybar_dat, divnorm_model, type="flex")
  foo <- rbind(compare_models(shrinky_fit, shrinky_fit, parind1 = 3),
               compare_models(div_count_fit, shrinky_fit),
               compare_models(div_flex_fit, shrinky_fit))
  foo <- cbind(response=response, model=c("Compressed", "Normalized: count", "Normalized: flex"), foo)
  foo[, pval:=formatC(pval, digits=1, format="g") |> str_replace("\\+0?", "\\+") |> str_replace("\\-0?", "\\-")]
  foo[, lr_pred:=formatC(lr_pred, digits=2, format="g") |> str_replace("\\+0?", "\\+") |> str_replace("\\-0?", "\\-")]
  return(foo)
}
