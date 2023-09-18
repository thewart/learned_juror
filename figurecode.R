library(rstanarm)
library(lme4)
library(cowplot)
library(tidybayes)
extrafont::loadfonts()
source("miscfunctions.R")
source("prep_stan.R")
evcond <- c("balanced","credible","defenseless")
ratecond <- c("without","with")

evidscheme <- c("Baseline"="black","Inculpatory"="#ba4040ff","Exculpatory"="#406bbaff", "Ambiguous"="#765884ff")
condscheme <- c("Balanced"="#994DA0","Credible"="#F87EBE", "Defenseless"="#E60026")
theme_set(theme_minimal_grid(font_size = 12, font_family = extrafont::choose_font("Arial")))
update_geom_defaults("pointrange",list(shape=19,size=0.5))
evidshape <- 16
lsz <- 14
figpath <- "../figures/"

load("../analysis/exp1_fitsplus.Rdata")
load("../analysis/exp2_fits.Rdata")
load("../analysis/exp2b_fits.Rdata")
load("../analysis/exp3_fits.Rdata")


#### Marginal impact of evidence ####
diff_plt <- ggplot(marginaldiff,aes(y=mean,x=evbal,ymin=lb,ymax=ub,color=type,shape=source,linetype=source)) + 
  geom_line(position = position_dodge(width=0.5)) + geom_pointrange(position = position_dodge(width=0.5),linetype=1) + 
  geom_hline(yintercept = 0) + ylab("Change in points") + xlab("Evidence balance") +
  scale_color_manual(values=evidscheme[c(2,3)], guide="none") + scale_shape_discrete("Source:") + scale_linetype_discrete("Source:") +
  coord_cartesian(clip = "off",ylim=c(-18,18),xlim=c(-3.1,3.1)) + 
  annotate("text",x=-2.5,y=-24,label="\u27F5 Exculpatory",fontface="italic",size=3) + 
  annotate("text",x=2.5,y=-24,label="Inculpatory \u27F6",fontface="italic",size=3)


flat <- data.table(evbal=c(-3,3,-3,3),mean=c(10,10,-10,-10),type=rep(c("Inculpatory","Exculpatory"),each=2), source="Behavior")
flat[,c("lb","ub"):=.(mean,mean)]
ggsave(paste0(figpath,"fltplt.pdf"), width = 8, height = 4, 
       delete_layers(diff_plt,"GeomPointrange") %+% flat + scale_linetype_discrete(guide="none") + 
         ylim(layer_scales(diff_plt)$y$get_limits()))

slope <- data.table(evbal=c(-3,3,-3,3),mean=c(5,15,-15,-5),type=rep(c("Inculpatory","Exculpatory"),each=2),source="Behavior")
slope[,c("lb","ub"):=.(mean,mean)]
ggsave(paste0(figpath,"slpplt.pdf"), width = 8, height = 4, 
       delete_layers(diff_plt,"GeomPointrange") %+% slope + scale_linetype_discrete(guide="none") +
         ylim(layer_scales(diff_plt)$y$get_limits()))

#### Case strength vs guilt probability ####
bvr <- stan_glmer(bardguilt ~ 1 + scale(rating) + (1 + scale(rating) | uid),iter=400,data=ldat,family = binomial)
foo <- apply(ldat[,((0:100)-mean(rating))/sd(rating)] %*%
               t(extract(bvr$stanfit,"scale(rating)")[[1]]), 1, function(x) x + extract(bvr$stanfit,"(Intercept)")[[1]])
foo <- post_summary_dt(1/(1+exp(-foo)))
foo$rating <- 0:100

bvr_c <- stan_glmer(bardguilt ~ 0 + cut_interval(rating,5) + (0 + cut_interval(rating,5) | uid),iter=1000,data=ldat,family = binomial)
juh <- 1/(1+exp(-extract(bvr_c$stanfit,"beta")[[1]])) %>% post_summary_dt()
juh$rating <- seq(10,90,20)
guiltvstrength_plt <- ggplot(foo, aes(y=mean,ymax=ub,ymin=lb,x=rating)) + geom_line() + geom_ribbon(alpha=0.25) + geom_pointrange(data=juh) +
  ylab("Guilty vote probability") + xlab("Case strength rating")


#### Evidence config PPC ####
bdat <- cbind(makebalancedat(ldat),ldat,Yhat=extract(rfit,"Yhat")[[1]] %>% colMeans())
pphat <- bdat[,.(Yhat=mean(Yhat),rating=mean(rating),evbal=unique(n_inculp-n_exculp)),
              by=.(physical,document,witness,character)]

evconfig_plt <- ggplot(pphat,aes(y=Yhat,x=rating)) + geom_point() + geom_abline() + 
  labs(x = "Predicted rating", y="Rating", fill= "Evidence balance")


# ggsave(paste0(figpath,"evconfig_bw.pdf"), width = 6, height = 4,
#        ggplot(pphat,aes(y=rating,x=Yhat)) + geom_point() + geom_abline() + 
#          labs(x = "Predicted rating", y="Rating") + geom_smooth(method="lm",color="black",linetype=2,se=F))

#### Evidence weights ####
effdt <- combine_ab(rfit,"mu_alpha_resp","mu_beta_resp") %>% post_summary_dt() %>% label_dt()
weights_plt <- ggplot(effdt,aes(y=mean,x=type,color=level)) + geom_pointrange(aes(ymin=lb,ymax=ub),position = position_dodge(width=0.3),shape=evidshape) +
  xlab("Evidence type") + ylab("Points") + geom_hline(yintercept = 0) +
  scale_color_manual("Evidence valence",breaks = names(evidscheme)[-1],values=evidscheme) + scale_x_discrete(drop=F)

##### scenario effects ####
samps <- extract(rfit,c("mu_alpha","mu_beta","alpha_scen","beta_scen","mu_scale","sigma_scale"))
sceneffdat <- rbind(with(samps,sweep(alpha_scen,1,mu_alpha,"+") |> apply(2,pointscale,mu_scale,sigma_scale)) |>
                      post_summary_dt(1:31) |> cbind(type="Baseline"),
                    with(samps,sweep(apply(beta_scen,c(1,2),mean),1,rowMeans(mu_beta),"+") |> apply(2,pointscale,mu_scale,sigma_scale,T)) |>
                      post_summary_dt(1:31) |> cbind(type="Evidence"))
setnames(sceneffdat,"level","scenario")

samps <- extract(rfit_cond$credible,c("mu_alpha","mu_beta","alpha_subj","beta_subj","mu_scale","sigma_scale"))
avebeta_subj <- with(samps,sweep(apply(beta_subj,c(1,2),mean),1,rowMeans(mu_beta),"+") |> apply(2,pointscale,mu_scale,sigma_scale,T))
alpha_subj <- with(samps,sweep(alpha_subj,1,mu_alpha,"+") |> apply(2,pointscale,mu_scale,sigma_scale))
cor.test(colMeans(avebeta_subj_resp),colMeans(alpha_subj_resp))
qplot(y=colMeans(avebeta_subj_resp),x=colMeans(alpha_subj_resp))
rowMeans(apply(samps$beta_subj,c(1,2),mean) * samps$alpha_subj) |> qplot(bins=100)

avebeta_dat <- post_summary_dt(avebeta_subj,1:ncol(avebeta_subj))
names(avebeta_dat) <- c("beta_hat","beta_lb","beta_ub","subj")
alpha_dat <- post_summary_dt(alpha_subj,1:ncol(alpha_subj))
names(alpha_dat) <- c("alpha_hat","alpha_lb","alpha_ub","subj")
alphabeta_dat <- alpha_dat[avebeta_dat,on="subj"]
ggplot(alphabeta_dat,aes(y=beta_hat,x=alpha_hat)) + geom_point() + geom_errorbar(aes(ymin=beta_lb,ymax=beta_ub),width=0,alpha=0.25) + geom_errorbarh(aes(xmin=alpha_lb,xmax=alpha_ub),height=0,alpha=0.25)


#### MEGAFIGURE ASSEMBLE ####
task_plt <- ggdraw() + draw_image(paste0(figpath,"task.png"))
statement_plt <- ggdraw() + draw_image(paste0(figpath,"missionstatement_abridged.png"))
# row1 <- plot_grid(task_plt,
#                   guiltvstrength_plt + theme(plot.margin = margin(1,1,1,1,unit = "lines")),
#                   evconfig_plt + theme(plot.margin = margin(1,1,1,1,unit="lines")),
#                   nrow=1,labels = "auto",label_size = lsz, rel_widths = c(1,.9,.9))
# row2 <- plot_grid(weights_plt + theme(legend.position = "top",plot.margin = margin(r=1,unit = "lines")),
#                   diff_plt + theme(legend.position = "top",plot.margin = margin(l=1,unit = "lines")),nrow = 1,labels = c("d","e"),label_size = lsz)

col1 <- plot_grid(statement_plt + theme(plot.margin = margin(1,0,1,0,unit = "lines")),
          task_plt + theme(plot.margin = margin(1,1,1,1,unit = "lines")),
          rel_heights=c(1,.75),ncol=1,labels="auto",label_size = lsz)

row1 <- plot_grid(guiltvstrength_plt + theme(plot.margin = margin(1,1,1,1,unit = "lines")),
                   evconfig_plt + theme(legend.position = "top",plot.margin = margin(1,1,1,1,unit = "lines")),
                  nrow=1,label_size = lsz, rel_widths = c(1,1))
col2 <- plot_grid(row1,weights_plt + theme(plot.margin = margin(r=1,l=1,unit = "lines")), labels=c("c","d"),ncol=1)

fig1 <- plot_grid(col1,col2,rel_widths = c(.5,1))

row1 <- plot_grid(statement_plt + theme(plot.margin = margin(0,1,0,1,unit = "lines")), 
                  task_plt + theme(plot.margin = margin(1,1,1,1,unit = "lines")),
                  guiltvstrength_plt + theme(plot.margin = margin(1,1,1,1,unit = "lines")),
                  nrow=1, labels = "auto")
row2 <- plot_grid(evconfig_plt + theme(legend.position = "top",plot.margin = margin(1,1,1,1,unit = "lines")),
                  weights_plt + theme(plot.margin = margin(r=1,l=1,unit = "lines")),
                  rel_widths = c(.5,1), nrow = 1, labels = c("d","e"))
fig1 <- plot_grid(row1,row2,nrow = 2)


ggsave(paste0(figpath,"figure1.pdf"),fig1,cairo_pdf,width = 9,height = 6)

##### Experiment 2 ####

#### Direct comparison ####

bdat_cond <- merge(ldat_cond,makebalancedat(ldat_cond),by=c("uid","scenario"))
bdat_cond[,nev:=n_exculp + n_inculp + n_ambig]
bdat_cond[,evconf:=paste0(physical,document,witness,character)]
head2head <- bdat_cond[,mean(rating),by=.(epoch=cut_interval(question,3),cond_evidence,evconf,nev)] |> dcast(epoch+evconf+nev ~ cond_evidence)
levels(head2head$epoch) <-  c("Early pre-instruction","Late pre-instruction","Post-instruction")
ggplot(head2head,aes(y=credible-balanced,x=nev)) + geom_point() + facet_wrap(vars(epoch),labeller = ) + geom_smooth(method="lm") + 
  xlab("Pieces of observed evidence") + ylab("Credible - balanced ratings (points)")


foo <- ldat_cond[,.(rating,scenario,cond_evidence,uid),.(paste0(physical,document,witness,character))]
juh <- foo[paste0 %in% foo[cond_evidence=="credible",unique(paste0)], mean(rating),by=.(paste0,cond_evidence)]
juh <- rbind(dcast(juh[cond_evidence!="defenseless"],paste0 ~ cond_evidence)[,type:="Credible"],
             dcast(juh[cond_evidence!="credible"],paste0 ~ cond_evidence)[!is.na(defenseless)][,type:="Defenseless"],use.names=F)
dcplt <- ggplot(juh,aes(y=credible-balanced,x=balanced,color=type)) + geom_point() + geom_hline(yintercept = 0) + geom_smooth(span=1,color="black") + 
  xlab("Balanced context (points)") + ylab("Unbalanced context (points)") + 
  scale_color_manual(NULL,values=condscheme[-1])

#### Time trend ####
bdat_cond <- merge(ldat_cond,makebalancedat(ldat_cond),by=c("uid","scenario"))

qbyq <- sapply(0:30, function(t) lm(rating ~ 0 + cond_evidence + cond_evidence:I(n_inculp-n_exculp),data=bdat_cond[question==t]) |> 
                 coef()) |> t() |> as.data.table()
qbyq$question <- 0:30
qint <- qbyq[,c(1:3,7),with=F]
qslo <- qbyq[,4:7,with=F]
setnames(qint,1:3,c("Balanced","Credible","Defenseless"))
setnames(qslo,1:3,c("Balanced","Credible","Defenseless"))
qint <- melt(qint,id.vars = "question",value.name = "est")
qslo <- melt(qslo,id.vars = "question",value.name = "est")
qslo$type <- "Evidence balance"
qint$type <- "Baseline"
qdat <- rbind(qint,qslo)

base <- matrix(0,nrow=12,ncol=31)
int_bal <- base
int_bal[1,] <- 1
int_bal[5,] <- 0:30
int_cred <- int_bal
int_cred[2,] <- 1
int_cred[8,] <- 0:30
int_def <- int_bal
int_def[3,] <- 1
int_def[9,] <- 0:30

slope_bal <- base
slope_bal[4,] <- 1
slope_bal[10,] <- 0:30
slope_cred <- slope_bal
slope_cred[6,] <- 1
slope_cred[11,] <- 0:30
slope_def <- slope_bal
slope_def[7,] <- 1
slope_def[12,] <- 0:30

mat2dat <- function(X,Y,func,condition=c("Balanced","Credible","Defenseless"),type,value.name) {
  dat <- lapply(X,func,Y) |> as.data.table()
  names(dat) <- condition
  dat$question <- 0:30
  dat <- melt(dat,id.vars="question",value.name = value.name)
  dat$type <- type
  return(dat)
}

trendat <- rbind(mat2dat(list(int_bal,int_cred,int_def),fixef(fit_trend),calc_con,type="Baseline",value.name="est"),
                 mat2dat(list(slope_bal,slope_cred,slope_def),fixef(fit_trend),calc_con,type="Evidence balance",value.name="est"))
sedat <- rbind(mat2dat(list(int_bal,int_cred,int_def),vcov(fit_trend),calc_con_se,type="Baseline",value.name="se"),
               mat2dat(list(slope_bal,slope_cred,slope_def),vcov(fit_trend),calc_con_se,type="Evidence balance",value.name="se"))

trendat <- trendat[sedat,on=.(question,variable,type)]
trendplt <- ggplot(trendat,aes(y=est,x=question,color=variable,fill=variable)) + geom_point(data=qdat,alpha=0.4) + geom_line() + 
  geom_ribbon(aes(ymin=est-se,ymax=est+se),alpha=0.4,linetype=0) + scale_color_brewer(NULL,palette = "Dark2") + 
  scale_fill_brewer(NULL,palette = "Dark2") + facet_wrap("type",scale="free_y") + 
  ylab("Points") + xlab("Case #")

#### Evidence weights by context ####
# e2maineff_samps <- (gather_draws(rfit_cond,mu_alpha_resp[cond],mu_beta_resp[cond,evidence]) |> 
#                       setDT() |> parse_evidence())[,cond:=set_factor_context(cond)] |> cull_defenseless()
# e2maineff <- rbind(weights_within_valence(e2maineff_samps),
#                    average_weights(e2maineff_samps)[cond!="Defenseless"])
# e2maineff[,notbaseline:=valence!="Baseline"]
# e2weights_plt <- weights_by_cond_plot(e2maineff,colorscale = condscheme) + facet_wrap(vars(notbaseline),ncol=1,scales="free_y",drop=F) + 
#   scale_y_continuous("Case strength (points)",breaks=seq(-10,70,10)) + theme(axis.text.x=element_text(angle=30,hjust=1),strip.text.x = element_blank(),panel.spacing = unit(1,"lines"))

e2subjeff <- subjeff_cond(rfit_cond,ldat_cond[,.(cond=cond_evidence[1]),by=uid])
e2maineff <- maineff_cond(rfit_cond)
ggplot(e2maineff,aes(y=mean,x=cond,color=cond)) + geom_point(data=e2subjeff,alpha=0.2,position = position_jitter(width = .1,height=0)) + geom_pointrange(aes(ymin=lb,ymax=ub,fill=cond),color="black",shape=22) + 
  facet_wrap(vars(valence),scales="free",nrow=1,strip.position = "bottom") + scale_x_discrete(NULL,breaks=NULL) + 
  scale_y_continuous("Evidence weight (points)",breaks=scales::extended_breaks(n=6)) + scale_color_manual("Condition",values=condscheme) + scale_fill_manual("Condition",values=condscheme)

#### Case strength distributions ####
foo <- merge(ldat_cond,makebalancedat(ldat_cond),by=c("uid","scenario"))
uc <- foo[,.(cond_evidence,physical,document,witness,character,evbal=n_inculp-n_exculp)] |> unique()
uc <- uc[,.N,by=.(cond_evidence,evbal)][,prop := N/sum(N),by=cond_evidence]
uc <- rbind(uc,data.table(cond_evidence=c("credible","defenseless"),evbal=c(-2,0),N=c(0,0),prop=c(0,0)))
ucplt <- ggplot(uc,aes(x=evbal,y=prop,color=str_to_sentence(cond_evidence))) + geom_line() + xlab("Evidence balance") + ylab("Proportion of cases") + scale_color_brewer(NULL,palette = "Dark2")

##### Megafigure 2! ####
fig2 <- plot_grid(ucplt + scale_colour_brewer(guide=NULL,palette = "Dark2") + theme(plot.margin = margin(1,2,0,1,"lines")),
                  e2weights_plt + theme(legend.position = "top"), 
                  dcplt + scale_color_manual(guide=NULL,values=RColorBrewer::brewer.pal(3,"Dark2")[-1]) + theme(plot.margin = margin(1,2,0,1,"lines")),
                  trendplt + scale_colour_brewer(guide=NULL,palette = "Dark2") + scale_fill_brewer(guide=NULL,palette = "Dark2"),
                  nrow=2,rel_widths = c(.6,1),align = "h",axis = "b",labels = c("a","c","b","d"),label_size = lsz)
ggsave(paste0(figpath,"figure2.pdf"),fig2,cairo_pdf,width = 9,height = 6)

#### Learning model ####
pltdat <- cbind(ldat_cond[,.(question,cond=cond_evidence)],
                Baseline=rlfit$par$baseline,
                Inculpatory=rlfit$par$inceff,
                Exculpatory=rlfit$par$exceff,
                Ambiguous=rlfit$par$ambeff,
                wbar=rlfit$par$Wbar,
                truncbonus=rlfit$par$trunc_bonus)
pltdat[,cond:=set_factor_context(cond)]

eldat <- pltdat[question %in% 0 | question %in% 30,lapply(.SD,mean),by=.(question>10,cond)] |> 
  melt(id.vars = c("question","cond"),variable.name = "valence",variable.factor=F)
eldat <- eldat[!(valence %in% c("wbar","truncbonus"))] |> cull_defenseless()
# setnames(eldat,"cond_evidence","cond")
# eldat <- aggsummlabels(eldat)
# eldat[(cond == "Defenseless" & level %in% c("Ambiguous","Exculpatory","Combined")),value:=NA]
eldat[,period:=(function(x) factor(c("Start","End"),levels = c("Start","End"))[x+1])(question)]

lcdat <- pltdat[,.(`Suppression penalty`=mean(pnorm(truncbonus/rlfit$par$scale)*100-50),
                   `Observed weight`=mean(pnorm(wbar/rlfit$par$scale)*100-50),
                   `Unobserved weight`=mean(pnorm((wbar+truncbonus)/rlfit$par$scale)*100-50)),
                by=.(question,cond)]
lcdat <- melt(lcdat,id.vars = c("question","cond"))
lcdat[,variable:=factor(variable,c("Unobserved weight","Observed weight","Suppression penalty"))]
# ggplot(lcdat,aes(y=value,x=question,color=str_to_title(cond))) + geom_line() + facet_wrap(vars(variable),nrow=1,scales="free_y") + 
#   scale_color_manual("Condition",breaks=c("Balanced","Credible","Defenseless"),values=condscheme) + xlab("Case #") + ylab("Case strength (points)")
ggplot(lcdat,aes(y=value,x=question,linetype=variable,color=cond)) + geom_line() + facet_wrap(vars(cond),nrow=1) + 
  scale_color_manual("Condition",values=condscheme) + xlab("Case #") + ylab("Case strength (points)") + geom_hline(yintercept = 0)

eldat[,value:=pnorm(value/rlfit$par$scale)*100]
eldat[valence!="Baseline",value:=value-50]
# eldat[,notbl:= level!="Baseline"]
eldat[,cond:=factor(cond,levels = c("Credible","Defenseless","Balanced"))]
eldat[,valence:=set_factor_valence(valence)]
pev <- ggplot(eldat[valence!="Baseline"],aes(y=value,color=cond,x=period,group=cond)) +geom_point(position = position_dodge(width=0.5),size = 2) + 
  geom_line(position = position_dodge(width=0.5)) + facet_wrap(vars(valence),strip.position = "bottom") + 
  ylab(NULL) + xlab(NULL) + scale_color_manual("Condition",breaks=c("Balanced","Credible","Defenseless"),values=condscheme)
pbase <- ggplot(eldat[valence=="Baseline"],aes(y=value,color=cond,x=period,group=cond)) + geom_point(position = position_dodge(width=0.5),size = 2) + 
  geom_line(position = position_dodge(width=0.5)) + facet_wrap(vars(valence),strip.position = "bottom") + 
  ylab("Effective weight (points)") + xlab(NULL) + scale_color_manual("Condition",breaks=c("Balanced","Credible","Defenseless"),values=condscheme,guide=NULL)

# prepostplt <- ggplot(eldat,aes(y=value,color=cond,x=period,shape=period,group=cond)) +
#   geom_point(position = position_dodge(width=0.5),size = 2) + geom_line(position = position_dodge(width=0.5)) + 
#   scale_color_manual(NULL,breaks=c("Balanced","Credible","Defenseless"),values = RColorBrewer::brewer.pal(3,"Dark2")) + 
#   scale_shape_manual(NULL,labels=c("Early","Late"), values=c(Start="square",End="triangle")) +
#   scale_y_continuous("Effective weight (points)") + facet_wrap(vars(level),scales = "free") + 
#   scale_x_discrete(NULL,label=NULL) + theme(strip.text.y = element_blank(),panel.spacing.x = unit(-1,"lines"), panel.spacing.y = unit(1,"lines"))
# 
# row1 <- plot_grid(ggdraw() + draw_image(paste0(figpath,"learner_juror.png")) + theme(plot.margin = margin(t=1,b=1,unit="lines")),
#                   ggdraw() + draw_image(paste0(figpath,"expectations.png")) + theme(plot.margin = margin(t=1,b=1,unit="lines")),
#                   rel_widths = c(.6,1), labels = "auto", label_size = lsz, nrow=1)
# row2 <- plot_grid(
#   lcplt + scale_colour_brewer(guide=NULL,palette = "Dark2") + theme(plot.margin = margin(t=2,l=2,r=2,unit = "lines")),
#   prepostplt + theme(legend.position = "right"),
#   labels = c("c","d"), nrow=1, rel_widths = c(0.6,1), label_size = lsz)

# fig3 <- plot_grid(row1,row2,ncol = 1)
  
# ggsave(paste0(figpath,"figure3.pdf"),fig3,cairo_pdf,width = 9,height = 6)

##### Binary choice exp 2a&b, 3 ####
foo <- gather_draws(bfit_rate,mu_alpha[cond],mu_beta[cond,evidence]) |> setDT() |> parse_evidence()
juh <- rbind(weights_within_valence(foo),average_weights(foo)) |> set_factor_rate() |> cull_defenseless()

moo <- (gather_draws(bfit_cond,mu_alpha[cond],mu_beta[cond,evidence]) |> setDT() |> 
          parse_evidence())[,cond:=set_factor_context(cond)] |> cull_defenseless()
buh <- rbind(weights_within_valence(moo),average_weights(moo)[cond!="Defenseless"])
buh$rating <- "Exp2"

lst <- list("Without" = "Exp. 2b: Without Rating",
            "With" = "Exp. 2b: With Rating",
            "Exp2" = "Exp. 2a")
explabeller <- function(variable,value) return(lst[value])
bard_plt <- weights_by_cond_plot(rbind(buh,juh),bard = T) + facet_wrap(vars(rating),ncol=1,labeller = explabeller)

ggsave("bardfig.pdf",bard_plt,width = 6,height = 6)
##### Exp 2b ratings ####
e3maineff_samps <- (gather_draws(rfit_rate,mu_alpha_resp[cond],mu_beta_resp[cond,evidence]) |> 
                      setDT() |> parse_evidence())[,cond:=set_factor_context(cond)] |> cull_defenseless()
e3maineff <- rbind(weights_within_valence(e3maineff_samps),
                   average_weights(e3maineff_samps)[cond!="Defenseless"])
e3maineff[,notbaseline:=valence!="Baseline"]
e3weights_plt <- weights_by_cond_plot(e3maineff) 
ggsave("2b_ratings.pdf",e3weights_plt,width = 6,height = 3)

#### Rating distributions ####
binsize <- 5
ppcdt <- ppc_ydist(ldat$rating,yrep,1)
ppcdt$y$type <- "Population"
ppcdt$yrep$type <- "Population"
popdist_plt <- ggplot(ppcdt$yrep,aes(x=rating,y=proportion/nrow(ldat))) + geom_ribbon(alpha=0.5,aes(ymin=lb/nrow(ldat),ymax=ub/nrow(ldat))) + 
  geom_line(linetype=2) + geom_col(data=ppcdt$y,alpha=0.5) + ylab("Proportion")  + facet_wrap("type") +
  scale_x_continuous(labels=c(0,25,50,75,100))

ppcinddt <- list(y=data.table(),yrep=data.table())
# ind <- sample.int(standat$Nsubj,9)
ind <- c(7,47,107)
for (i in ind) {
  iind <- standat$Subj==i
  ippcdt <- ppc_ydist(ldat$rating[iind],yrep[,iind],binsize)
  ippcdt$y$Subj <- paste0("S",i)
  ippcdt$yrep$Subj <- paste0("S",i)
  ppcinddt$y <- rbind(ppcinddt$y,ippcdt$y)
  ppcinddt$yrep <- rbind(ppcinddt$yrep,ippcdt$yrep)
}
inddist_plt <- ggplot(ppcinddt$yrep,aes(x=rating,y=proportion/31)) + geom_line(linetype=2) + 
  geom_col(data=ppcinddt$y,alpha=0.5) + ylab("Frequency") + facet_wrap("Subj",ncol=1,strip.position = "right") + 
  scale_x_continuous(labels=c(0,25,50,75,100)) + ylab(NULL) + scale_y_continuous(name = NULL)

ppc_plt <- plot_grid(popdist_plt + xlab(NULL) + theme(plot.margin = margin(c(5.5,5.5,20,5.5))),
                     inddist_plt + xlab(NULL) + ylab(NULL) + theme(plot.margin = margin(c(5.5,5.5,20,5.5)))) + 
  draw_label("Rating", y=0.03, size = 11)

ggsave("rppc.pdf",ppc_plt,width=6,height = 4)

#### Pairwise interaction vs linear model comparison ####
standat <- makestandat(ldat, interact = T)
intvsmain <- function(standat,fitint) {
  intconts <- str_split_fixed(colnames(standat$Z),":",n=2) %>% apply(1,function(y) colnames(standat$X) %in% y)
  foo <- extract(fitint,"mu_beta")[[1]] %*% intconts %>% post_summary_dt()
  colnames(foo) <- c("compsum","clb","cub")
  return(cbind(foo,extract(fitint,"mu_lambda")[[1]] %>% post_summary_dt))
}
intdt <- intvsmain(standat,rfitint)

intconts <- str_split_fixed(colnames(standat$Z),":",n=2) %>% apply(1,function(y) colnames(standat$X) %in% y)
foo <- extract(rfitint,"mu_beta")[[1]] %*% intconts
juh <- extract(rfitint,"mu_lambda")[[1]]
buh <- cbind(post_summary_dt(foo),post_summary_dt(juh+foo))
colnames(buh) <- colnames(intdt)
nonlinear_plt <- ggplot(buh,aes(x=compsum,y=mean)) + geom_point() + geom_errorbar(aes(ymin=lb,ymax=ub),alpha=0.25,width=0) +
  geom_errorbarh(aes(xmin=clb,xmax=cub),alpha=0.25,height=0) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + 
  ylab("Total effect (a.u.)") + xlab("Sum of main effects (a.u.)") + geom_abline() 

ggsave("nonlinear.pdf",nonlinear_plt,width = 6,height=4)


##### Capstone evidence conditions ####

# Rating
cap_weights_samps <- gather_draws(rfit_cap,mu_alpha_pre_resp[cond],mu_beta_pre_resp[cond,evidence]) |> setDT() |> parse_evidence()
cap_weights_samps[,cond:=set_factor_context(cond)]

cap_cond_contrasts <- (cap_weights_samps[,mean(.value),by=.(cond,valence,.draw)] |> 
                         dcast(valence + .draw ~ cond, value.var="V1"))[,post_summary_dt(Balanced-Credible),by=valence]

cap_weights_summary <- rbind(weights_within_valence(cap_weights_samps),average_weights(cap_weights_samps))
cap_weights_plt <- weights_by_cond_plot(cap_weights_summary)

# BARD
cap_weights_bard_samps <- gather_draws(bfit_cap,mu_alpha[cond],mu_beta[cond,evidence]) |> setDT() |> parse_evidence()
cap_weights_bard_samps[,cond:=set_factor_context(cond)]

cap_weights_bard_summary <- rbind(weights_within_valence(cap_weights_bard_samps),average_weights(cap_weights_bard_samps))
cap_weights_bard_plt <- weights_by_cond_plot(cap_weights_bard_summary,bard = T)


cap_cond_bcontrasts <- (cap_weights_bard_samps[,mean(.value),by=.(cond,valence,.draw)] |> 
                         dcast(valence + .draw ~ cond, value.var="V1"))[,post_summary_dt(Balanced-Credible),by=valence]

cap_weights_bard_plt <- weights_by_cond_plot(cap_weights_bard_summary,T)

##### capstone intervention effects ####
# Rating

e3subjeff <- rbind(
  cbind(subjeff_cond(rfit_cap,ldat_cap[,.(cond=cond_evidence[1]),by=uid],average_over="all",post_int=T), capped=T),
  cbind(subjeff_cond(rfit_cap,ldat_cap[,.(cond=cond_evidence[1]),by=uid],average_over="all"), capped=F)
)

e3maineff <- rbind(
  cbind(maineff_cond(rfit_cap,post_int=T,average_over = "all"), capped=T),
  cbind(maineff_cond(rfit_cap,post_int=F,average_over = "all"), capped=F)
)

ggplot(e3maineff,aes(y=mean,x=capped,color=cond)) + geom_line(aes(group=cond),color="black",position = position_dodge(width=.5)) + 
  geom_pointrange(aes(ymin=lb,ymax=ub,fill=cond),color="black",shape=22,position = position_dodge(width=.5)) +
  geom_point(data=e3subjeff,alpha=0.2,position = position_jitterdodge(dodge.width = .5)) +
  facet_wrap(vars(valence),scales="free",strip.position = "bottom") + scale_color_manual("Condition",values = condscheme[-3]) + scale_fill_manual("Condition",values = condscheme[-3])


# cap_intervene_samps <- (gather_draws(rfit_cap,mu_alpha_pre_resp[cond],mu_alpha_post_resp[cond],mu_beta_ave_pre_resp[cond],mu_beta_ave_post_resp[cond]) |> 
#                        setDT())[,cond:=set_factor_context(cond)]
# cap_intervene_summary <- cap_intervene_samps[,post_summary_dt(.value),by=.(cond,.variable)] |> parse_capstone()
# cap_intervene_summary[,type:=factor(type,levels = unique(type))]
# capcrossovers_plt <- ggplot(cap_intervene_summary,aes(y=mean,ymin=lb,ymax=ub,x=capped,color=cond,group=cond)) + geom_pointrange(position = position_dodge(width=.5)) + 
#   facet_wrap(vars(type),scales="free",strip.position = "bottom") + geom_line(position = position_dodge(width=.5)) + scale_color_brewer(NULL,palette = "Dark2") +
#   scale_x_discrete(NULL,labels=c("Pre-inst.","Post-inst.")) + ylab("Case strength (points)") + theme(strip.placement = "outside")

# BARD
lambda_samps <- gather_draws(bfit_cap,mu_lambda[cond,type]) |> setDT()
lambda_samps[,cond:=set_factor_context(cond)]
lambda_samps[type==1,valence:="Baseline"]
lambda_samps[type==2,valence:="Ave. evidence"]

cap_intervene_bard_samps <- rbind(cap_weights_bard_samps[valence=="Baseline",-4,with=F],average_weights(cap_weights_bard_samps,F))
cap_intervene_bard_samps <- cap_intervene_bard_samps[lambda_samps[,.(cond,.draw,valence,capval=.value)],on=.(cond,valence,.draw)]
cap_intervene_bard_summary <- rbind(
  cbind(cap_intervene_bard_samps[,post_summary_dt(.value),by=.(cond,valence)],capped=F),
  cbind(cap_intervene_bard_samps[,post_summary_dt(.value+capval),by=.(cond,valence)],capped=T))

ggplot(cap_intervene_bard_summary,aes(y=mean,ymin=lb,ymax=ub,x=capped,color=cond,group=cond)) + geom_pointrange(position = position_dodge(width=.5)) + 
  facet_wrap(vars(valence),scales="free",strip.position = "bottom") + geom_line(position = position_dodge(width=.5)) + scale_color_brewer(NULL,palette = "Dark2") +
  scale_x_discrete(NULL,labels=c("Pre-inst.","Post-inst.")) + ylab("Guilt judgment (log odds)") + theme(strip.placement = "outside")

##### Direct comparison v2 ####
bdat_cap <- merge(ldat_cap,makebalancedat(ldat_cap),by=c("uid","scenario"))
bdat_cap[,nev:=n_exculp + n_inculp + n_ambig]
bdat_cap[,evconf:=paste0(physical,document,witness,character)]
head2head <- (bdat_cap[,mean(rating),by=.(epoch=cut_interval(question,3),cond_capstone,cond_evidence,evconf,nev)] |> 
                dcast(epoch+cond_capstone+evconf+nev ~ cond_evidence))[!is.na(credible)]
levels(head2head$epoch) <-  c("Early pre-instruction","Late pre-instruction","Post-instruction")
ggplot(head2head,aes(y=credible-balanced,x=nev)) + geom_point() + facet_wrap(vars(epoch),labeller = ) + geom_smooth(method="lm") + 
  xlab("Pieces of observed evidence") + ylab("Credible - balanced ratings (points)")


##### Mini-mega figure ####
fig4 <- plot_grid(ggdraw() + draw_image(paste0(figpath,"task3.png")) + theme(plot.margin = margin(1,1,1,1,unit="lines")),
          capcrossovers_plt + theme(legend.position = "top"),
          nrow=1, rel_widths = c(.75,1), labels = "auto", label_size = lsz)

ggsave(paste0(figpath,"figure4.pdf"),fig4,cairo_pdf,width = 9,height = 3)
