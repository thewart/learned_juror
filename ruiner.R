whichtask <- "exculpatory_conditional"
source("process_data.R")
ldat_cond <- makelegaldat(scendat,subjdat,clickdat)
setkey(ldat_cond,uid,question)

ldat_cond[,evconf:=paste0(physical,document,witness,character)]
ldat_cond_cred <- ldat_cond[evconf %in% ldat_cond[cond_evidence=="credible",evconf]]

standat <- makestandat(ldat_cond_cred,cond="cond_evidence",levelref=T)
rfit_cond_cred <- sampling(condmodel,standat,chains=4,iter=500,warmup=200,pars=c(basepars,rpars))
rfit_cond_cred <- recover_types(rfit_cond_cred,standat$ref)

ldat_cond[,evconf:=paste0(physical,document,witness,character)]
ldat_cond_def <- ldat_cond[evconf %in% ldat_cond[cond_evidence=="defenseless",evconf]]

standat <- makestandat(ldat_cond_def,cond="cond_evidence",levelref=T)
rfit_cond_def <- sampling(condmodel,standat,chains=4,iter=500,warmup=200,pars=c(basepars,rpars))
rfit_cond_def <- recover_types(rfit_cond_def,standat$ref)

foo <- rbind(
  cbind(maineff_cond(rfit_cond),fit_on="All"),
  cbind(maineff_cond(rfit_cond_cred),fit_on="Credible only"),
  cbind(maineff_cond(rfit_cond_def)[valence %in% c("Baseline","Inculpatory")],fit_on="Defenseless only")
)
