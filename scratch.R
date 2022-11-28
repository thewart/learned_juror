whichtask <- "exculpatory_conditional"
source("~/code/casereveal/analysis/process_data.R")
source("~/code/casereveal/analysis/prep_stan.R")
source("~/code/casereveal/analysis/funcscratch.R")
ldat <- makelegaldat(scendat,subjdat,clickdat)
setkey(ldat,uid,question)
ordcond <- c("random","credible","inculpatory")

bdat <- makebalancedat(ldat) %>% merge(ldat,by=c("uid","scenario"))
setkey(bdat,uid,question)

mvplt <- rbind(mvwnft(bdat[cond_evidence=="random"],5,"random"),
               mvwnft(bdat[cond_evidence=="credible"],5,"credible"),
               mvwnft(bdat[cond_evidence=="inculpatory"],5,"inculpatory"))
# meltplt <- melt(mvplt,id.vars=c("question","condition"),measure.vars=c("bl","ev","exceff"),value.name = "points")
# meltplt$se <- melt(mvplt,id.vars=c("question","condition"),measure.vars=c("bl_se","inceff_se","exceff_se"))$value
meltplt <- melt(mvplt,id.vars=c("question","condition"),measure.vars=patterns("[lf]$"),value.name = "points")
meltplt$se <- melt(mvplt,id.vars=c("question","condition"),measure.vars=patterns("se"))$value

ggplot(meltplt,aes(y=points,x=question,color=condition,fill=condition)) + geom_line() + geom_ribbon(aes(ymin=points-se,ymax=points+se),alpha=0.25,linetype=0) +
  facet_wrap("variable",scales = "free")



standat <- makestanerdat(ldat)
model <- stan_model("/home/seth/code/legalmodels/refcase_transfunc1_flat.stan")
standat$Cond[standat$Cond==3] <- 1
guh <- optimizing(model,standat,as_vector=F)

pltdat <- cbind(ldat[,.(question,cond_evidence)],baseline=guh$par$baseline,inculpatory=guh$par$inceff,
                exculpatory=guh$par$exceff,ambiguous=guh$par$ambeff,combined=guh$par$combined)
eldat <- pltdat[question %in% range(question),lapply(.SD,mean),by=.(question,cond_evidence)] %>% 
  melt(id.vars = c("question","cond_evidence"),variable.name = "level")
setnames(eldat,"cond_evidence","evidence")
eldat <- aggsummlabels(eldat)
eldat[(evidence == "inculpatory" & level %in% c("ambiguous","exculpatory","combined")),value:=0]
eldat[,period:=(function(x) factor(c("Start","End"),levels = c("Start","End"))[x+1])(question/30)]
ggplot(eldat,aes(x=level,y=value,fill=evidence)) + geom_col(position = position_dodge()) + facet_wrap("period")


ggplot(pltdat,aes(y=exceff,x=question,color=cond_evidence)) + geom_line()