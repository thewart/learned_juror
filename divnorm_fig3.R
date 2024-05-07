##### binary ####
ybar_bin <- ybarify(bdat, resp="binary")

shrinky_max_onesided <- fit_onesided(ybar_bin, shrinky_beta)
plot_ypred_poly(ybar_bin, shrinky_max_onesided$par$ypred, 3)
lm_est_int(ybar_bin, shrinky_max_onesided$par$ypred, 3) |> summary()
test_int(ybar_bin, shrinky_max_onesided$par$ypred, 3)

shrinky_max <- fit_full(ybar_bin, shrinky_beta)
plot_ypred_poly(ybar_bin, shrinky_max$par$ypred, 1)
lm_est_int(ybar_bin, shrinky_max$par$ypred, 1) |> summary()
test_int(ybar_bin, shrinky_max$par$ypred, 1)

divglm <- fit_divglm(ybar_bin, shrinky_beta_div, type="flex")
kdiff <- purrr::list_c(divglm$par[1:5]) |> length() - purrr:::list_c(shrinky_max$par[1:3]) |> length()
lldiff <- divglm$value - shrinky_max$value
1-pchisq(2*lldiff, kdiff)
exp(lldiff - kdiff) 
plot_ypred_poly(ybar_bin, divglm$par$yhat, 1)

##### rating ####
ybar_dat <- ybarify(bdat, resp="rating")

shrinky_max_onesided <- fit_onesided(ybar_dat, shrinky)
p1 <- plot_ypred_poly(ybar_dat, shrinky_max_onesided$par$ypred, degree = 3)
lm_est_int(ybar_dat, shrinky_max_onesided$par$ypred, 3) |> summary()
lm_est_resid(ybar_dat, shrinky_max_onesided$par$ypred, 3) |> summary()
test_int(ybar_dat, shrinky_max_onesided$par$ypred, 3)

shrinky_max <- fit_full(ybar_dat, shrinky)
p2 <- plot_ypred_poly(ybar_dat, shrinky_max$par$ypred, 1)
lm_est_int(ybar_dat, shrinky_max$par$ypred, 1) |> summary()
test_int(ybar_dat, shrinky_max$par$ypred, 1)

##### funky variance sanity check ####
juh <- bdat[,.(rating, nev, evconf)][
  data.table(evconf=ybar_dat$evconf, ybar=ybar_dat$ybar, ypred=shrinky_max$par$ypred), on="evconf"][
    , .(se=sd(rating/100)), by=.(nev, ypred, ybar, evconf)]
# ggplot(juh[nev>0], aes(y=se, x=ypred, color=factor(nev))) + geom_point() + geom_smooth(method="lm", formula=y ~ poly(x,2), se=F)
juh[, ypred_bin:=cut(ybar, breaks=c(0, .25, .4, .6, 1))]
# juh[, ypred_bin:=cut_number(ypred,4)]
juh[, cor.test(se, nev), by=.(ypred_bin)]
pb1 <- ggplot(juh[nev>0], aes(y=se, x=nev, color=ordered(nev))) + geom_smooth(aes(color=NULL), method="lm", se=F, color="black") + geom_jitter(width=0.1, height=0) + 
  facet_wrap(vars(ypred_bin), ncol=2, strip.position = "top") + ylab("Rating SD") + xlab("") + scale_color_viridis_d(end=0.9)
anova(lm(ybar-ypred ~ scale(se) + poly(ypred,1) + poly(ypred,1):scale(se), data=juh[nev>0]),
      lm(ybar-ypred ~ scale(se) + nev + poly(ypred,1) + poly(ypred,1):scale(se) + poly(ypred,1):nev, data=juh[nev>0]))

divglm <- fit_divglm(ybar_dat, shrinky_div, type="linear_count")
p3 <- plot_ypred_poly(ybar_dat, divglm$par$yhat, 1)
test_int(ybar_dat, divglm$par$yhat, 1)

kdiff <- purrr::list_c(divglm$par[1:5]) |> length() - purrr:::list_c(shrinky_max$par[1:3]) |> length()
lldiff <- divglm$value - shrinky_max$value
1-pchisq(2*lldiff, kdiff)
exp(lldiff - kdiff) # AIC

##### illustrative plotz
p4 <- ggplot(NULL) + geom_function(fun=pdiv, args=list(pars=divglm$par, nev=1), aes(color="1"), linewidth=0.75) + 
  geom_function(fun=pdiv, args=list(pars=divglm$par, nev=2), aes(color="2"), linewidth=0.75) + 
  geom_function(fun=pdiv, args=list(pars=divglm$par, nev=3), aes(color="3"), linewidth=0.75) + 
  geom_function(fun=pdiv, args=list(pars=divglm$par, nev=4), aes(color="4"), linewidth=0.75) + 
  scale_color_viridis_d("Evidence \n count", end=0.9) + ylab("Predicted case strength") + 
  scale_x_continuous("Sum of weights", limits = c(-4, 4), labels=c("", "\U2190 Exculpatory", "", "Inculpatory \U2192", ""))
# pt <- plot_grid(p1 + theme(legend.position = "none") + xlab("Compressed sum of weights \n (single-valence fit)"), 
#           p2 + scale_y_continuous(NULL, labels=NULL) + theme(legend.position = "none") + xlab("Compressed sum of weights \n (all cases fit)"), 
#           p3 + scale_y_continuous(NULL, labels=NULL) + xlab("Normed sum of weights \n (all cases fit)"),
#           align="h", axis="tblr", nrow=1, rel_widths=c(1.05, 1, 1.25), labels = c("a", "", ""))
pt <- plot_grid(p1 + theme(legend.position = "none") + xlab("Summed-compressed predictions \n (single-valence fit)"), 
                p2 + scale_y_continuous(NULL, labels=NULL) + xlab("Summed-compressed predictions \n (all cases fit)"), 
                p3 + theme(legend.position = "none") + xlab("Summed-normed-compressed predictions \n (all cases fit)"),
                p4,
                align="h", axis="tblr", rel_widths = c(0.85, 1), nrow=2, labels="auto")

pb <- plot_grid(pb1 + theme(legend.position = "none"), p4, rel_widths = c(1.1, 1), labels=c("b","c"))
plot_grid(pt,pb, ncol=1, rel_heights = c(1.2, 1))

alligator <- bdat[(n_exculp==0 | n_inculp==0) & (character=="none"), .(ybar=mean(rating/100), se=sd(rating/100)/sqrt(.N)), by=.(n_inculp, n_exculp, n_ambig)]
alligator <- rbind(alligator[n_inculp==1][, `:=` (valence="Inculpatory", n_clear=n_inculp)],
                   alligator[n_exculp==1][, `:=` (valence="Exculpatory", n_clear=n_exculp)]) |> cbind(source="Observed")
# p1 <- nomnom_alligator(alligator[(n_inculp<2) & (n_exculp<2)])
lmer(rating/100 ~ n_ambig*n_inculp + n_ambig*n_exculp + (1|uid), 
     data=bdat[(character=="none") & ((n_exculp + n_inculp)==1)]) |> summary()

# alligator_cond <- bdat_cond[(cond_evidence %in% c("balanced", "credible")) & (n_exculp==0 | n_inculp==0) & (character=="none"),
#                             .(ybar=mean(rating/100), se=sd(rating/100)/sqrt(.N)), by=.(n_inculp, n_exculp, n_ambig)]
# nomnom_alligator(alligator_cond[(n_inculp<2) & (n_exculp<2)])
# lmer(rating/100 ~ n_ambig*n_inculp + n_ambig*n_exculp + (1|uid), 
#      data=bdat_cond[(cond_evidence %in% c("balanced", "credible")) & (character=="none") & ((n_exculp + n_inculp) %in% c(0,1))]) |> summary()

ybar_dat$yhat <- shrinky_max$par$yhat
alligator_shrinky <- ybar_dat[(n_exculp==0 | n_inculp==0) & (character=="none"), .(ybar=mean(yhat), se=0), by=.(n_inculp, n_exculp, n_ambig)]
alligator_shrinky <- rbind(alligator_shrinky[n_exculp==1][, `:=` (valence="Exculpatory", n_clear=n_exculp)],
                   alligator_shrinky[n_inculp==1][, `:=` (valence="Inculpatory", n_clear=n_inculp)]) |> cbind(source="Summation-compression")
# p2 <- nomnom_alligator(alligator_shrinky[(n_inculp<2) & (n_exculp<2)])

divglm <- fit_divglm(ybar_dat, shrinky_div, type="flex")
ybar_dat$yhat <- divglm$par$yhat
alligator_div <- ybar_dat[(n_exculp==0 | n_inculp==0) & (character=="none"), .(ybar=mean(yhat), se=0), by=.(n_inculp, n_exculp, n_ambig)]
alligator_div <- rbind(alligator_div[n_exculp==1][, `:=` (valence="Exculpatory", n_clear=n_exculp)],
                           alligator_div[n_inculp==1][, `:=` (valence="Inculpatory", n_clear=n_inculp)]) |> cbind(source="Sum-normalize-compress")
# p3 <- nomnom_alligator(alligator_div[(n_inculp<2) & (n_exculp<2)])
alligator <- rbind(alligator, alligator_shrinky, alligator_div)

ggplot(alligator, aes(y=ybar, x=n_ambig, ymin=ybar-2*se, ymax=ybar+2*se, color=valence)) + 
  facet_grid(factor(valence, levels=unique(valence)) ~ factor(source, levels=unique(source)), scales = "free_y") + 
  geom_smooth(method="lm", se=F, color="black") + geom_pointrange() + theme(legend.position = "none") +
  scale_color_manual("Evidence valence",breaks = names(evidscheme)[-1],values=evidscheme) + scale_x_continuous("Ambiguoius evidence count", breaks=seq(0,2)) + scale_y_continuous("Case strength rating", breaks=scales::breaks_extended(n=4))

# nomnom_alligator <- function(dt) {
#   alligator <- rbind(dt[n_exculp==0][, `:=` (valence="Inculpatory", n_clear=n_inculp)],
#                      dt[n_inculp==0][, `:=` (valence="Exculpatory", n_clear=n_exculp)])
#   plt <- ggplot(alligator[valence=="Inculpatory"], aes(y=ybar, x=n_clear, color=ordered(n_ambig))) + 
#     geom_pointrange(aes(ymin=ybar-2*se, ymax=ybar+2*se), position = position_dodge(width=0.2)) + geom_line(position = position_dodge(width=0.2)) + 
#     geom_pointrange(data=alligator[valence=="Exculpatory"], aes(ymin=ybar-2*se, ymax=ybar+2*se), position = position_dodge(width=0.2)) + 
#     geom_line(data=alligator[valence=="Exculpatory"], aes(y=ybar, x=n_clear, color=ordered(n_ambig)), position = position_dodge(width=0.2), linetype=2) + 
#     xlab("# clear evidence") + ylab("Case strength rating") + scale_color_ordinal("# ambig. \n evidence")
#   return(plt)
# }
# 
# plot_grid(p1 + theme(legend.position = "none") + ylim(c(0.1,0.75)) + scale_x_continuous(breaks=c(0,1)) + ggtitle("Observed") + theme(plot.title = element_text(size = 14, face="plain")),
#           p2 + theme(legend.position = "none") + ylab(NULL) + ylim(c(0.1,0.75)) + scale_x_continuous(breaks=c(0,1)) + ggtitle("Summation-compression") + theme(plot.title = element_text(size = 14, face="plain")),
#           p3 + ylab(NULL) + ylim(c(0.1,0.75)) + scale_x_continuous(breaks=c(0,1)) + ggtitle("Sum-normalize-compress") + theme(plot.title = element_text(size = 14, face="plain")), 
#           nrow=1, rel_widths=c(1,1,1.2))

# logistic_beta <- stan_model("models/logistic_beta.stan")
# ybar_singleton <- ybar_dat[nev %in% c(0,1)]
# X <- model.matrix(~ physical + document + witness + character, ybar_singleton)[,-1]
# X_pred <- model.matrix(~ physical + document + witness + character, ybar_dat)[,-1]
# standat <- list(N=nrow(X), P=ncol(X), X=X, ybar=ybar_singleton$ybar, M=ybar_singleton$M,
#                 X_pred=X_pred, N_pred=nrow(X_pred))
# max_singleton <- optimizing(logistic_beta, standat, as_vector=F, verbose=T)
# 
# ybar_comp <- cbind(ybar_dat, yhat=max_singleton$par$ypred)
# ybar_comp[,lm(ybar-yhat ~ yhat*I(n_inculp+n_exculp+n_ambig))] |> summary()
# cross <- ybar_comp[nev==0, ybar]
# abdat <- ybar_comp[nev>0, .(b=lm(ybar-yhat ~ 0 + I(yhat-cross))$coef), by=nev][, a:= -cross*b]
# ggplot(ybar_comp[nev>1], aes(y=ybar, x=yhat, color=factor(nev))) + geom_point(alpha=0.3) +
#   geom_abline(data=abdat[nev>1], aes(intercept=a, slope=b+1, color=factor(nev)), linewidth=0.75) + geom_abline() + xlab("Sum of evidence weights (compressed)") + ylab("Observed guilt prob.") +
#   scale_color_discrete("Evidence count")


