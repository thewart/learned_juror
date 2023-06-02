library(rstanarm)
library(lme4)
library(cowplot)
source("miscfunctions.R")
source("prep_stan.R")
evcond <- c("balanced","credible","defenseless")
ratecond <- c("without","with")

load("../analysis/exp1_fitsplus.Rdata")
load("../analysis/exp2_fits.Rdata")
load("../analysis/exp3_fits.Rdata")

colorscheme <- c("Baseline"="black","Inculpatory"="red","Exculpatory"="blue", "Ambiguous"="magenta")
theme_set(theme_minimal_grid(font_size = 10, font_family = "arial"))
lsz <- 11
figpath <- "../analysis/presentation_figs/"

#### Marginal impact of evidence ####
diff_plt <- ggplot(marginaldiff,aes(y=mean,x=evbal,ymin=lb,ymax=ub,color=type,shape=source,linetype=source)) + 
  geom_line(position = position_dodge(width=0.5)) + geom_pointrange(position = position_dodge(width=0.5),linetype=1) + 
  geom_hline(yintercept = 0) + ylab("Change in points") + xlab("Evidence balance") +
  scale_color_manual(values=colorscheme[c(2,3)], guide="none") + scale_shape_discrete("Source:") + scale_linetype_discrete("Source:") +
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


ggsave(paste0(figpath,"evconfig_bw.pdf"), width = 6, height = 4,
       ggplot(pphat,aes(y=rating,x=Yhat)) + geom_point() + geom_abline() + 
         labs(x = "Predicted rating", y="Rating") + geom_smooth(method="lm",color="black",linetype=2,se=F))

#### Evidence weights ####
effdt <- combine_ab(rfit,"mu_alpha_resp","mu_beta_resp") %>% post_summary_dt() %>% label_dt()
weights_plt <- ggplot(effdt,aes(y=mean,x=type,color=level)) + geom_pointrange(aes(ymin=lb,ymax=ub),position = position_dodge(width=0.3)) +
  xlab("Evidence type") + ylab("Points") + geom_hline(yintercept = 0) +
  scale_color_manual("Evidence level:",breaks = names(colorscheme)[-1],values=colorscheme) + scale_x_discrete(drop=F)

# weights_plt %+% effdt[level=="inculpatory"] + ylim(layer_scales(weights_plt)$y$get_limits())
# ggsave(paste0(figpath,"subweights1.pdf"), width = 8, height = 4,
#        weights_plt %+% effdt[level=="inculpatory"] + ylim(layer_scales(weights_plt)$y$get_limits()))
# ggsave(paste0(figpath,"subweights2.pdf"), width = 8, height = 4,
#        weights_plt %+% effdt[level %in% c("inculpatory","exculpatory")] + ylim(layer_scales(weights_plt)$y$get_limits()))
# ggsave(paste0(figpath,"subweights3.pdf"), width = 8, height = 4,
#        weights_plt %+% effdt[level!="baseline"] + ylim(layer_scales(weights_plt)$y$get_limits()))
# ggsave(paste0(figpath,"weights.pdf"), width = 8, height = 4, weights_plt)


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
task_plt <- ggdraw() + draw_image("task.png")
row1 <- plot_grid(task_plt,
                  guiltvstrength_plt + theme(plot.margin = margin(1,1,1,1,unit = "lines")),
                  evconfig_plt + theme(plot.margin = margin(1,1,1,1,unit="lines")),
                  nrow=1,labels = "auto",label_size = lsz, rel_widths = c(1,.9,.9))
row2 <- plot_grid(weights_plt + theme(legend.position = "top",plot.margin = margin(r=1,unit = "lines")),
                  diff_plt + theme(legend.position = "top",plot.margin = margin(l=1,unit = "lines")),nrow = 1,labels = c("d","e"),label_size = lsz)

fig1 <-plot_grid(row1,row2,ncol = 1)

ggsave("figure1.pdf",fig1,cairo_pdf,width = 9,height = 6)

#### Direct comparison ####
foo <- ldat_cond[,.(rating,scenario,cond_evidence,uid),.(paste0(physical,document,witness,character))]
juh <- foo[paste0 %in% foo[cond_evidence=="credible",unique(paste0)], mean(rating),by=.(paste0,cond_evidence)]
juh <- rbind(dcast(juh[cond_evidence!="defenseless"],paste0 ~ cond_evidence)[,type:="Credible"],
             dcast(juh[cond_evidence!="credible"],paste0 ~ cond_evidence)[!is.na(defenseless)][,type:="Defenseless"],use.names=F)
dcplt <- ggplot(juh,aes(y=credible,x=balanced,color=type)) + geom_point() + geom_abline() + 
  xlab("Balanced context (points)") + ylab("Unbalanced context (points)") + 
  scale_color_manual(NULL,values=RColorBrewer::brewer.pal(3,"Dark2")[-1])

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
e2maineff <- aggsumm_cond(rfit_cond,evcond,resp_scale = T)[!(condition == "defenseless" & level %in% c("ambiguous","exculpatory","combined"))]
# e2maineff[level=="baseline",c("mean","lb","ub"):=.(mean-50,lb-50,ub-50)]
e2maineff <- aggsummlabels(e2maineff)

e2maineff[,notbaseline:=level!="Baseline"]
e2weights_plt <- ggplot(e2maineff, aes(x=level,y=mean,color=condition)) +
  geom_pointrange(aes(ymin=lb,ymax=ub),position = position_dodge(width=0.3)) + geom_hline(data=data.table(y=c(0,NA),notbaseline=c(T,F)),aes(yintercept=y)) +
  xlab(NULL) + scale_color_brewer(NULL,palette = "Dark2") + scale_y_continuous("Case strength (points)",breaks=seq(-10,70,10)) + 
  facet_wrap(vars(notbaseline),ncol=1,scales="free_y",drop=F) +  
  theme(axis.text.x=element_text(angle=30,hjust=1),strip.text.x = element_blank(),panel.spacing = unit(1,"lines")) 
  # e2maineff[,notbaseline:=NULL]

# ggsave(paste0(figpath,"condweights.pdf"), width = 8, height = 4, condweights)

# ggsave(paste0(figpath,"condweights_baseline.pdf"), width = 8, height = 4,
# condweights %+% e2maineff[level=="baseline"] + scale_x_discrete(drop=F) + ylim(layer_scales(condweights)$y$get_limits()))

e2maineff_samps <- (gather_draws(rfit_cond_alt,mu_alpha_resp[cond],mu_beta_resp[cond,evidence]) |> 
                      setDT() |> parse_evidence())[,cond:=set_factor_context(cond)] |> cull_defenseless()
e2maineff <- e2maineff_samps[,mean(.value),by=.(cond,valence,.draw)][,post_summary_dt(V1),by=.(cond,valence)]
e2maineff[,notbaseline:=valence!="Baseline"]

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
ggsave("figure2.pdf",fig2,cairo_pdf,width = 9,height = 6)

#### Learning model ####
pltdat <- cbind(ldat_cond[,.(question,condition=cond_evidence)],Baseline=rlfit$par$baseline,Inculpatory=rlfit$par$inceff,
                Exculpatory=rlfit$par$exceff,Ambiguous=rlfit$par$ambeff,Combined=rlfit$par$combined,wbar=rlfit$par$Wbar)

eldat <- pltdat[question %in% range(question),lapply(.SD,mean),by=.(question,condition)] |> 
  melt(id.vars = c("question","condition"),variable.name = "level",variable.factor=F)
eldat <- eldat[level != "wbar"]
# setnames(eldat,"cond_evidence","condition")
eldat <- aggsummlabels(eldat)
eldat[(condition == "Defenseless" & level %in% c("Ambiguous","Exculpatory","Combined")),value:=NA]
eldat[,period:=(function(x) factor(c("Start","End"),levels = c("Start","End"))[x+1])(question/30)]

lcdat <- pltdat[,.(wbar=mean(pnorm(wbar/rlfit$par$scale)*100-50)),by=.(question,condition)]
lcplt <- ggplot(lcdat,aes(y=wbar,x=question,color=str_to_title(condition))) + geom_line() + 
  scale_color_brewer(NULL,palette = "Dark2") + xlab("Case #") + ylab("Expected strength (points)")


eldat[,value:=pnorm(value/rlfit$par$scale)*100]
eldat[level!="Baseline",value:=value-50]
eldat[,notbl:= level!="Baseline"]
eldat[,condition:=factor(condition,levels = c("Credible","Defenseless","Balanced"))]
prepostplt <- ggplot(eldat[(level!="Combined")],aes(y=value,color=condition,x=period,shape=period,group=condition)) +
  geom_point(position = position_dodge(width=0.5),size = 2) + geom_line(position = position_dodge(width=0.5)) +
  geom_hline(data=data.table(y=c(0,NA),notbl=c(T,F)),aes(yintercept=y)) + 
  scale_color_manual(NULL,breaks=c("Balanced","Credible","Defenseless"),values = RColorBrewer::brewer.pal(3,"Dark2")) + 
  scale_shape_manual(NULL,labels=c("Learning Start","Learing End"), values=c(Start="square",End="triangle")) +
  scale_y_continuous("Case strength (points)",breaks=seq(-10,70,10)) + facet_grid(notbl ~ level,switch = "x",scales = "free_y", space = "free") + 
  scale_x_discrete(NULL,label=NULL) + theme(strip.text.y = element_blank(),panel.spacing.x = unit(-1,"lines"), panel.spacing.y = unit(1,"lines"))

fig3 <- plot_grid(ggdraw() + draw_image(magick::image_read_svg("lineup_alt.svg")),
                    ggdraw() + draw_image(magick::image_read_svg("boxtheory.svg")),
                    lcplt + scale_colour_brewer(guide=NULL,palette = "Dark2") + theme(plot.margin = margin(t=2,l=2,r=2,unit = "lines")),
                    nrow = 1, rel_widths = c(1,0.6,1), labels = "auto", label_size = lsz) |> 
          plot_grid(prepostplt + theme(legend.position = "right"),labels = c("","d"), nrow = 2,scale = c(1,.9),label_size = lsz)
ggsave("figure3.pdf",fig3,cairo_pdf,width = 9,height = 6)

##### Binary choice exp 2a&b, 3 ####
postsumm <- rbind(aggsumm_cond(bfit_rate[[ratecond[1]]],evcond,"condition","rating","Exp. 2b: Without Rating"),
                  aggsumm_cond(bfit_rate[[ratecond[2]]],evcond,"condition","rating","Exp. 2b: With Rating"),
                  aggsumm_cond(bfit_cond,evcond,"condition","rating","Exp. 2a")
                  )[!(condition == "defenseless" & level %in% c("ambiguous","exculpatory","combined"))]
postsumm <- aggsummlabels(postsumm)
# postsumm <- postsumm[!(level=="Combined" & condition=="Defenseless")]
# postsumm <- postsumm[(level %in% c("Baseline","Combined")) | (level=="Inculpatory")]
# postsumm[level == "Combined",level:="all evidence"]
bard_plt <- ggplot(postsumm,aes(x=level,y=mean,ymin=lb,ymax=ub,color=condition)) + geom_pointrange(position = position_dodge(width=.25)) +
  facet_wrap("rating",ncol=1) + xlab(NULL) + ylab("Guilty judgement (log odds)") + scale_color_discrete(NULL) +
  geom_hline(yintercept=0) + theme(axis.text.x=element_text(angle=30,hjust=1),legend.position = "top")
ggsave("bardfig.pdf",bard_plt,width = 6,height = 6)
##### Exp 2b ratings ####
e3maineff <- aggsumm_cond(rfit_rate,evcond,resp_scale = T)[!(condition == "defenseless" & level %in% c("ambiguous","exculpatory","combined"))]
# e2maineff[level=="baseline",c("mean","lb","ub"):=.(mean-50,lb-50,ub-50)]
e3maineff <- aggsummlabels(e3maineff)

e3maineff[,notbaseline:=level!="Baseline"]
e3weights_plt <- ggplot(e3maineff, aes(x=level,y=mean,color=condition)) +
  geom_pointrange(aes(ymin=lb,ymax=ub),position = position_dodge(width=0.3)) + geom_hline(data=data.table(y=c(0,NA),notbaseline=c(T,F)),aes(yintercept=y)) +
  xlab(NULL) + scale_color_brewer(NULL,palette = "Dark2") + scale_y_continuous("Case strength (points)") + 
  theme(axis.text.x=element_text(angle=30,hjust=1),legend.position = "top") 
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
