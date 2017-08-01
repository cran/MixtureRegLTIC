###############################
### function: mixtureLogitAFT()
###############################
MixtureLogitAFT=function(formula,eventprobreg=~1,locationreg=~1,scalereg=~1,var.entry,var.mixturetype=NULL,var.weight=NULL,
                data,time.origin=0,shape=NULL){

  cat("Running ...")

  time1=Sys.time()

  call=match.call()

  ##################
  ##### Parameterize
  ##################
  shapereg=~1
  out=parameterize(formula,eventprobreg,locationreg,scalereg,shapereg,data,shape,var.entry)
  dataset=out$dataset
  covariates=out$covariates
  var.centime=out$var.centime
  var.truntime=out$var.truntime
  is.reg=out$is.reg
  x=out$x
  fix.par=out$fix.par
  par=out$par
  is.interval=out$is.interval
  is.truncation=out$is.truncation
  minmax.time=out$minmax.time
  var.obsweight=var.weight
  is.weight=F; if(!is.null(var.obsweight)) is.weight=T

  ###################
  ##### Error message
  ###################
  if(is.null(dataset)) return("Error: data !!")
  if(is.null(var.centime)) return("Error: survival time !!")
  if(any(is.na(match(var.centime,names(dataset))))) return("Error: survival time !!")
  if(!is.null(var.truntime)){
    if(any(is.na(match(var.truntime,names(dataset))))) return("Error: truncation time !!")
  }
  if(!is.null(covariates$names[-1])){
    if(any(is.na(match(covariates$names[-1],names(dataset))))) return("Error: covariates !!")
  }
  if(!is.null(var.mixturetype)){
    if(any(is.na(match(var.mixturetype,names(dataset))))) return("Error: variable of mixture type !!")
  }
  if(!is.null(var.obsweight)){
    if(any(is.na(match(var.obsweight,names(dataset))))) return("Error: variable of observation weight !!")
  }

  ###################
  ##### Keep variable
  ###################
  dataset=dataset[,c(var.centime,var.truntime,covariates$names[-1],var.mixturetype,var.obsweight)]

  #########################
  ##### Keep observations
  #########################
  dataset=obsKeep(dataset,covariates,var.centime,var.truntime)

  num.obs=nrow(dataset)

  ### The information of groups stratified by covariates
  ### "groups.obs" (befor add interaction)
  if (length(covariates$names[-1])>0){
    groups.obs=groupObs(dataset[,covariates$names[-1]],covariates$levels)
  }else{
    groups.obs=list()
    groups.obs$labels="0";groups.obs$levels=0;groups.obs$levelsTrue=T;groups.obs$obs=rep(0,dim(dataset)[1])
  }
  if(length(covariates$labels)==length(groups.obs$labels))groups.obs$labels=covariates$labels
  groups.obs$total.obs=numeric(length(groups.obs$levels)) ### "groups.obs$total.obs"
  groups.obs$case.obs=numeric(length(groups.obs$levels)) ### "groups.obs$case.obs"
  buffer.dataset=cbind(dataset[,var.centime[3]],groups.obs$obs)
  total.obs=table(buffer.dataset[,2])
  groups.obs$total.obs[match(as.numeric(names(total.obs)),groups.obs$levels)]=total.obs
  index=which(buffer.dataset[,1] != 2)
  if(length(index)>0) buffer.dataset=buffer.dataset[index,]
  case.obs=table(buffer.dataset[,2])
  groups.obs$case.obs[match(as.numeric(names(case.obs)),groups.obs$levels)]=case.obs

  ## Create "levels.mixturetype" from "var.mixturetype"
  levels.mixturetype=NULL
  if(!is.null(var.mixturetype)){
    buffer=unique(cbind(dataset[,var.mixturetype],groups.obs$obs))
    buffer=buffer[order(buffer[,2]),]
    if(is.null(dim(buffer))){
      buffer=buffer[!is.na(buffer[2])]
      levels.mixturetype=buffer[1]
    }else{
      buffer=buffer[!is.na(buffer[,2]),]
      levels.mixturetype=buffer[,1]
    }
  }

  ############
  ### "weight"
  ############
  weight=rep(1,num.obs)
  if(!is.null(var.obsweight)) weight=dataset[,var.obsweight]

  ##############################################
  ### "mtype" from "var.mixturetype" or "is.reg"
  ##############################################
  if(!is.null(var.mixturetype)){
    mtype=dataset[,var.mixturetype]
  }else{
    if(is.reg$beta==F){
      mtype=rep(1,num.obs)       # 1: one component
    }else{
      mtype=rep(2,num.obs)       # 2: two components with cure
    }
  }

  ##############
  ### "survtime"
  ##############
  ### var.centime (1:exact; 2:right censor; 3:interval censor; 4:left censor)
  ### var.truntime (1:no truncation; 2:left truncation; 3:interval truncation; 4:right truncation)
  var=var.centime
  if(is.truncation) var=c(var,var.truntime)
  survtime=dataset[,var]
  index=which(survtime[,var.centime[3]]==2)
  if(length(index)>0) survtime[index,var.centime[2]]=NA
  index=which(survtime[,var.centime[3]]==4)
  if(length(index)>0) survtime[index,var.centime[1]]=NA
  if(is.truncation){
    index=which(survtime[,var.truntime[3]]==2)
    if(length(index)>0) survtime[index,var.truntime[2]]=NA
    index=which(survtime[,var.truntime[3]]==4)
    if(length(index)>0) survtime[index,var.truntime[1]]=NA
  }

  #########################################
  ### Adjust survival time from time origin
  #########################################
  var=c(var.centime[1:2])
  if(is.truncation) var=c(var,var.truntime[1:2])
  survtime[,var]=survtime[,var]-time.origin

  ######################################
  ### Take nature log with survival time
  ######################################
  survtime[,var]=log(survtime[,var])

  ###############################
  ##### "fix.par" => "is.fix.par"
  ###############################
  is.fix.par=list()
  if (length(covariates$q) != length(fix.par$q)){
    fix.par$q=NULL
  }
  is.fix.par$q = !is.null(fix.par$q)

  ############################
  ##### Initial Value => "par"
  ############################
  survtime.turnbull=cbind(survtime,weight=weight)
  index=which(survtime.turnbull[,3]==1)
  if(length(index)>0) survtime.turnbull[index,2]=survtime.turnbull[index,1]
  if (is.truncation){
    turnbullEst=turnbullR(dataset=survtime.turnbull,var.centime=names(survtime.turnbull)[1:3],var.truntime=names(survtime.turnbull)[4:6],
             var.obsweight=names(survtime.turnbull)[7],is.graph=F,tolerance=1e-3)
    TurnbullEST=condition.turnbullR(turnbullEst,is.graph=F)
  }else{
    turnbullEst=turnbullR(dataset=survtime.turnbull,var.centime=names(survtime.turnbull)[1:3],var.truntime=NULL,
             var.obsweight=names(survtime.turnbull)[4],is.graph=F,tolerance=1e-3)
    TurnbullEST=condition.turnbullR(turnbullEst,is.graph=F)
  }
  TurnbullEST$p0=1-turnbullEst[1,"survival"]
  TurnbullEST$p2=turnbullEst[nrow(turnbullEst)-1,"survival"]
  if (TurnbullEST$SkewY > 0) TurnbullEST$SkewY=-1 else TurnbullEST$SkewY=1
  y=list(E=NULL,Var=NULL,Skew=NULL)
  y[1:3]=TurnbullEST[2:4]
  #
  if (is.reg$beta){
    par$beta[1]=log((1-TurnbullEST$p0-TurnbullEST$p2)/TurnbullEST$p2)
  }
  w=list(E=NULL,Var=NULL)
  w$E=(log(y$Skew^2)+digamma(y$Skew^{-2}))/y$Skew
  w$Var=trigamma(y$Skew^{-2})/y$Skew^2
  par$alpha[1]=log(sqrt(y$Var/w$Var))
  par$gamma[1]=y$E-w$E*exp(par$alpha[1])
  par$q[1]=y$Skew
  #
  if (is.fix.par$q) par$q=fix.par$q

  ##############
  ### "par.nfix"
  ##############
  par.nfix=NULL
  if (is.reg$beta){
    par.nfix=c(par.nfix,par$beta)
  }
  par.nfix=c(par.nfix,par$gamma,par$alpha)
  if (!is.fix.par$q){
    par.nfix=c(par.nfix,par$q)
  }

  #######################
  ### parameter estimator
  #######################
  options(warn=-1)
  convergence=1
  times.call=1; times.stop=3
  while(convergence!=0 & times.call<=times.stop){
    convergence=0
    parscale=rep(1^(times.call-1),length(par.nfix))

    cat(" ...")

    est=optim(par=par.nfix,fn=LLF,gr=GLLF,method="BFGS",control=list(parscale=parscale,reltol=1e-16),hessian=F,
        par.all=par,is.fix.par=is.fix.par,survtime=survtime,x=x,is.reg=is.reg,weight=weight,mtype=mtype)

    par.nfix=est$par

    ######################################
    ### statistic of regression result (1)
    ######################################
    gradient=-GLLF(par.nfix,par,is.fix.par,survtime,x,is.reg,weight,mtype)
    hessian=HLLF(par.nfix,par,is.fix.par,survtime,x,is.reg,weight,mtype)
    COV=NULL
    std=rep(NA,length(par.nfix))
    if(!any(is.na(hessian) | abs(hessian)==Inf)){
      if(det(hessian)>0){
        COV=solve(hessian,tol=1e-25)
        std=sqrt(diag(COV))
      }
    }
    if(est$convergence!=0){
      convergence=est$convergence
    }else if(any(is.na(hessian) | abs(hessian)==Inf)){
      convergence=2
    }else if(det(hessian)<=0){
      convergence=3
    }else if(any(is.na(std))){
      convergence=4
    }else if(max(abs(gradient))>1e-03){
      convergence=5
    }
    times.call=times.call+1
  }
  cat("\n")
  if(convergence==0){
    cat("Convergence"); cat("\n")
  }else{
    cat("Not convergence"); cat("\n")
  }

  ######################################
  ### statistic of regression result (2)
  ######################################
  chisq=rep(NA,length(par.nfix))
  pvalue=rep(NA,length(par.nfix))
  if(all(!is.na(std))){
    chisq=(par.nfix/std)^2
    pvalue=1-pchisq(chisq,1)
  }
  LogLF=-est$value
  AIC=2.0*(-LogLF+length(par.nfix))

  #######################
  ### "par.nfix" => "par"
  #######################
  start=1
  if (is.reg$beta){
    len=length(par$beta); par$beta=par.nfix[seq(start,start+len-1)]; start=start+len
  }
  len=length(par$gamma); par$gamma=par.nfix[seq(start,start+len-1)]; start=start+len
  len=length(par$alpha); par$alpha=par.nfix[seq(start,start+len-1)]; start=start+len
  if (!is.fix.par$q){
    len=length(par$q); par$q=par.nfix[seq(start,start+len-1)]; start=start+len
  }

  ############
  ### residual
  ############
  residuals.data=residualData(par,survtime,x)
  index=which(residuals.data[,3]==1)
  if(length(index)>0) residuals.data[index,2]=residuals.data[index,1]
  if(is.truncation){
    index=which(residuals.data[,4]==-Inf)
    if(length(index)>0){
      residuals.data[index,4]=NA
      residuals.data[index,6]=0
    }
  }
  vars=c(covariates$main,var.obsweight)
  if(length(vars)>0){
    res.dataset=cbind(residuals.data,dataset[vars])
  }else{
    res.dataset=residuals.data
  }

  #################################
  # right censored subjects: 2 => 0
  #  left censored subjects: 4 => 2
  #################################
  index=which(res.dataset[,3]==2)  # 2 => 0 (right censored subjects)
  if(length(index)>0) res.dataset[index,3]=0
  index=which(res.dataset[,3]==4)  # 4 => 2 (left censored subjects)
  if(length(index)>0) res.dataset[index,3]=2

  groups.obs=groups.obs[c("labels","total.obs","case.obs")]

  time2=Sys.time()
  run.time=time2-time1
  run.time=paste(as.numeric(run.time),units(run.time),sep=" ")

  cat(run.time); cat("\n")

  MixtureRegEST=list(call=call,convergence=convergence,covariates=covariates,time.origin=time.origin,is.truncation=is.truncation,minmax.time=minmax.time,
                     is.interval=is.interval,is.weight=is.weight,fix.par=fix.par,par=par,parest=par.nfix,
                     std=std,chisq=chisq,pvalue=pvalue,gradient=gradient,LLF=LogLF,AIC=AIC,COV=COV,groups.obs=groups.obs,mixturetype=levels.mixturetype,
                     res.dataset=res.dataset,run.time=run.time,method="mixtureLogitAFT")

  class(MixtureRegEST) = "mixture"

  return(MixtureRegEST)
}
if(F){
  par.all=par
  print(LLF(par.nfix,par.all,is.fix.par,survtime,x,is.reg,weight,mtype),16)
  GLLF(par.nfix,par.all,is.fix.par,survtime,x,is.reg,weight,mtype)
  HLLF(par.nfix,par.all,is.fix.par,survtime,x,is.reg,weight,mtype)
}
