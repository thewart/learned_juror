source("miscfunctions.R")
source("prep_stan.R")
source("divnorm_helpers.R")
library(lme4)
library(kableExtra)
library(cowplot)
theme_set(cowplot::theme_minimal_grid())
evidscheme <- c("Baseline"="black","Inculpatory"="#ba4040ff","Exculpatory"="#406bbaff", "Ambiguous"="#765884ff")
theme_set(cowplot::theme_minimal_grid(font_size=12))

shrinky_beta <- stan_model("models/shrinkage_t_beta.stan")
shrinky_beta_div <- stan_model("models/shrinkaget_beta_divnorm_logscale.stan")

shrinky <- stan_model("models/shrinkage_t_maxlik.stan")
shrinky_div <- stan_model("models/shrinkaget_divnorm_logscale.stan")

bdat <- readindat("exculpatory_burdenofproof", pthresh=0.05)
bdat_cond <- readindat("exculpatory_conditional", pthresh=0.05)[cond_evidence %in% c("balanced", "credible")]
bdat_rate <- readindat("exculpatory_rateless", pthresh=0.05)[cond_evidence %in% c("balanced", "credible")]

t1 <- rbind(maketable(ybarify(bdat, resp="rating"), shrinky, shrinky_div, "Case strength"),
      maketable(ybarify(bdat, resp="binary"), shrinky_beta, shrinky_beta_div, "Guilt judgment"))

##### test for equivalence ####
bdat_cond_credconf <- bdat_cond[evconf %in% bdat_cond[cond_evidence=="credible", unique(evconf)]]
anova(lmer(rating/100 ~ cond_evidence + evconf + (0 + cond_evidence|evconf) + (1|uid), data=bdat_cond_credconf, REML=F),
      lmer(rating/100 ~ evconf + (1|uid), data=bdat_cond_credconf, REML=F))
t2 <- rbind(maketable(ybarify(bdat_cond, resp="rating"), shrinky, shrinky_div, "Case strength"),
      maketable(ybarify(bdat_cond, resp="binary"), shrinky_beta, shrinky_beta_div, "Guilt judgment"))

good <- bdat_rate[,cor(n_inculp-n_exculp, as.numeric(bardguilt)), by=uid][V1>0 & !is.na(V1), uid]
bdat_rate <- bdat_rate[uid %in% good]
t3 <- rbind(maketable(ybarify(bdat_rate[cond_rating=="with"], resp="rating"), shrinky, shrinky_div, "Case strength"),
            maketable(ybarify(bdat_rate[cond_rating=="with"], resp="binary"), shrinky_beta, shrinky_beta_div, "Guilt w/ rating"),
            maketable(ybarify(bdat_rate[cond_rating=="without"], resp="binary"), shrinky_beta, shrinky_beta_div, "Guilt w/o rating"),
            maketable(ybarify(bdat_rate, resp="binary"), shrinky_beta, shrinky_beta_div, "Guilt (all)"))

bdat_cap <- readindat("capstone_prolific")[cond_evidence %in% c("balanced", "credible") & cond_capstone==F & (n_inculp+n_exculp+n_ambig)>0] 
capbar <- ybarify(bdat_cap, resp="rating")
capbar <- capbar[!is.na(se)]
rbind(maketable(capbar, shrinky, shrinky_div, "Case strength"),
      maketable(ybarify(bdat_cap, resp="binary"), shrinky_beta, shrinky_beta_div, "Guilt judgment"))

rbind(t1, t2, t3) |> kable("latex", digits=2, booktabs=T) |> collapse_rows(1, latex_hline = "major") |> 
  pack_rows("Experiment 1", 1, 6) |> pack_rows("Experiment 2", 7, 12) |> pack_rows("Experiment 3", 13, 21)


foo <- data.table(x=1:4, y=c(summary(lm(rating/100 ~ evconf, data=bdat))$r.squared,
                             summary(lm(rating/100 ~ uid, data=bdat))$r.squared,
                             summary(lm(rating/100 ~ scenario, data=bdat))$r.squared,
                             summary(lm(rating/100 ~ evconf + uid + scenario, data=bdat))$r.squared))


bvr <- stan_glmer(bardguilt ~ 1 + scale(rating) + (1 + scale(rating) | uid),iter=400,data=bdat,family = binomial)
foo <- apply(bdat[,((0:100)-mean(rating))/sd(rating)] %*%
               t(extract(bvr$stanfit,"scale(rating)")[[1]]), 1, function(x) x + extract(bvr$stanfit,"(Intercept)")[[1]])
foo <- post_summary_dt(1/(1+exp(-foo)))
foo$rating <- 0:100

bvr_c <- stan_glmer(bardguilt ~ 0 + cut_interval(rating,7) + (0 + cut_interval(rating,7) | uid),iter=1000,data=bdat, family = binomial)
juh <- 1/(1+exp(-extract(bvr_c$stanfit,"beta")[[1]])) %>% post_summary_dt()
juh$rating <- seq(0, 100, length.out=15)[seq(2, 15, 2)]
rvb_plt <- ggplot(foo, aes(y=mean,ymax=ub,ymin=lb,x=rating)) + geom_line() + geom_ribbon(alpha=0.25) + geom_pointrange(data=juh) +
  ylab("Guilt probability") + scale_x_continuous("Case strength rating (single response)", labels=c(0, 0.25, 0.5, 0.75, 1))

rvb_case_plt <- ggplot(NULL, aes(y=ybarify(bdat, resp = "binary")$ybar, x=ybarify(bdat, resp = "rating")$ybar)) + 
  geom_point(alpha=0.3) + geom_smooth(method="lm", formula=y ~ poly(x,1), color="black", se=F) + 
  xlab("Case strength rating (case average)") + ylab("Guilt probability") + ylim(c(0,1)) + xlim(c(0,1))

juh <- data.table(x=1:4, y=c(summary(lm(rating/100 ~ evconf, data=bdat))$r.squared,
summary(lm(rating/100 ~ uid, data=bdat))$r.squared,
summary(lm(rating/100 ~ scenario, data=bdat))$r.squared,
summary(lm(rating/100 ~ evconf + uid + scenario, data=bdat))$r.squared))
r2comp <- ggplot(juh, aes(y=y, x=x)) + geom_bar(stat = "identity") + 
  ylab(expression("R"^"2")) + scale_x_continuous(NULL, breaks=1:4, labels = c("Evidence set", "Participant", "Crime", "Combined"), minor_breaks = NULL)

