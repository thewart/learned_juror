makelegaldat <- function(scendat,subjdat,clickdat,fthresh=NULL) {
  msubj <- subjdat[!str_detect(uid,"test") & task_complete,uid]
  msubj <- msubj[msubj %in% subjdat[,.N,by=uid][N==1,uid]] #exclude subjects who repeated the task
  msubj <- msubj[msubj %in% clickdat[,.N==31,by=uid][V1==T,uid]] #exclude subjects with wrong number of reponses
  scendat <- scendat[uid %in% msubj]
  # clickdat <- clickdat[!is.na(guilty)]
  
  if ("legalguilt" %in% names(clickdat)) names(clickdat)[which(names(clickdat) == "legalguilt")] <- "bardguilt"
  
  dat <- merge(scendat,clickdat,by=c("uid","scenario"))
  
  # were there conditions?
  cond <- str_subset(colnames(subjdat),"cond")
  if (length(cond)>0) dat <- merge(dat,subjdat[,c("uid",cond),with=F],by="uid")
  # relable evidence conditions 
  if ("cond_evidence" %in% cond) {
    dat[cond_evidence == "random", cond_evidence:="balanced"]
    dat[cond_evidence == "inculpatory", cond_evidence:="defenseless"]
  }
  if ("cond_rating" %in% cond) {
    dat$cond_rating <- as.character(dat$cond_rating)
    dat[cond_rating == "0", cond_rating:="without"]
    dat[cond_rating == "1", cond_rating:="with"]
    dat[cond_rating == "without", rating:=NA]
  }
  # are there trial types?
  wcond <- str_subset(colnames(scendat),"cond")
  # if (length(wcond)>0) dat <- merge(dat,scendat[,c("uid","scenario",wcond),with=F],by=c("uid","scenario"))
  if ("cond_capstone" %in% wcond) dat[,cond_capstone:=(cond_capstone==1)]
  
  # orderize the evidence levels 
  evidord <- c("none", "clear_ex", "ambiguous", "clear_in")
  dat[,physical:=factor(physical,evidord)]
  dat[,document:=factor(document,evidord)]
  dat[,witness:=factor(witness,evidord)]
  dat[,character:=factor(character,evidord[-3])]
  # dat[, severity := scenario/max(scenario)]
  
  dat[,start:=as.POSIXct(start)]
  dat[,stop:=as.POSIXct(stop)]
  
  if (!is.null(fthresh)) {
    dat <- merge(dat,dat[,.(fst=anova(glm(bardguilt ~ rating,family=binomial))[2,2]),by=uid])
    dat <- dat[fst>fthresh]
  }
  return(dat)
}

makestandat <- function(dat,interact=F,binresp=F) { 
  Nsubj <- length(unique(dat$uid))
  Nscen <- length(unique(dat$scenario))
  
  # get upper and lower-bounded censored data
  L <- min(dat$rating)
  U <- max(dat$rating)
  
  # get censoring data frame
  X <- model.matrix(~ 1 + physical + document + witness + character, data=dat)[,-1]
  S <- dat[,factor(uid) %>% as.numeric]
  C <- dat$scenario
  
  # useful dimensions
  N <- dim(X)[1]
  P <- dim(X)[2]
  
  standat <- list( Nsubj=Nsubj, Nscen=Nscen, N=N, P=P, Scen=C, Subj=S, X=X)
  
  if (interact) {
    Z <- model.matrix(~ 1 + physical + document + witness + character + 
                        physical:document + physical:witness + physical:character + 
                        document:witness + document:character + witness:character, data=dat)[,-c(1:(ncol(X)+1))]
    Z <- Z[,which(apply(Z,2,uniqueN)==2)]
    standat <- c(standat, Z = list(Z), P2 = ncol(Z))
  }
  if (binresp) {
    standat <- c(standat,Y = list(dat[,as.numeric(bardguilt)]))
  } else {
    Y <- dat$rating
    # cens <- (Y == U) - (Y == L)
    # standat <- c(standat, L=L, U=U, Y=list(Y), cens=list(cens), D=1)
    standat <- c(standat, L=L, U=U, Y=list(Y), D=1)
  }
  
  return(standat)
}

makebalancedat <- function(dat) {
  return(dat[,.(n_exculp=str_detect(sapply(.(physical,document,witness,character),as.character),"ex") %>% sum(),
                                         n_inculp=str_detect(sapply(.(physical,document,witness,character),as.character),"in") %>% sum(),
                                         n_ambig=str_detect(sapply(.(physical,document,witness,character),as.character),"amb") %>% sum()),
                                         by=.(uid,scenario)])
}

makebalancedat_nochar <- function(dat) {
  return(dat[,.(n_exculp=str_detect(sapply(.(physical,document,witness),as.character),"ex") %>% sum(),
                n_inculp=str_detect(sapply(.(physical,document,witness),as.character),"in") %>% sum(),
                n_ambig=str_detect(sapply(.(physical,document,witness),as.character),"amb") %>% sum()),
             by=.(uid,scenario)])
}
