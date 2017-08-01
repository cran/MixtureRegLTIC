############################
### function: plotMixture()
############################
plotMixture=function(fit,dist="overall",curve="survival",
  xlab=NULL,ylab=NULL,main=NULL,col=NULL,lty=NULL,lwd=1,axes=T){

  est=c(fit,curvesLogAFT(fit))

  plotNPMLEsurv(est=est,dist=dist,curve=curve,type="l",
  xlab=xlab,ylab=ylab,main=main,col=col,lty=lty,lwd=lwd,axes=axes)
}

############################
### function: curvesLogAFT()
############################
curvesLogAFT=function(fit){

  n.group=length(fit$groups.obs$labels)

  buffer=apply(as.matrix(fit$groups.obs$labels),1,strsplit,split=",")
  matrix.cov=matrix(1,length(buffer),1+length(buffer[[1]][[1]]))
  for(k in 1:length(buffer)) matrix.cov[k,]=c(1,as.numeric(buffer[[k]][[1]]))

  ### "p1", "xbeta", "xgamma", "xalpha", "xq"
  p1=as.matrix(rep(1,dim(matrix.cov)[1]))
  if (!is.null(fit$covariates$beta)){
    xbeta=as.matrix(matrix.cov[,fit$covariates$beta])%*%as.vector(fit$par$beta)
    p1=exp(xbeta)/(1+exp(xbeta))
    if(!is.null(fit$mixturetype)){
      index=which(fit$mixturetype==1)
      if(length(index)>0) p1[index,]=1
    }
  }
  xgamma=as.matrix(matrix.cov[,fit$covariates$gamma])%*%as.vector(fit$par$gamma)
  xalpha=as.matrix(matrix.cov[,fit$covariates$alpha])%*%as.vector(fit$par$alpha)
  if(is.null(fit$fix.par$q)){
    xq=as.matrix(matrix.cov[,fit$covariates$q])%*%as.vector(fit$par$q)
  }else{
    xq=matrix(fit$fix.par$q,n.group)
  }

  ### "surv"
  time=seq(fit$minmax.time[1],fit$minmax.time[2],by=0.01)
  y=log(time-fit$time.origin)

  surv=list()
  for (i in 1:n.group){
    w=(y-xgamma[i])/exp(xalpha[i])
    q=rep(xq[i],length(w))

    surv1=as.data.frame(matrix(NA,length(time),7))
    names(surv1)=c("time","survival","survivalD","density","densityD","hazard","hazardD")

    surv1[,"time"]=time
    surv1[,"densityD"]=dw.fun(w,q)/(time*exp(xalpha[i]))
    surv1[,"survivalD"]=Sw.fun(w,q)
    surv1[,"hazardD"]=surv1[,"densityD"]/surv1[,"survivalD"]
    index=which(is.na(surv1[,"hazardD"]))
    if(length(index)>0){
      if(min(index)>1) surv1[index,"hazardD"]=surv1[min(index)-1,"hazardD"]
    }

    surv1[,"density"]=p1[i]*surv1[,"densityD"]
    surv1[,"survival"]=p1[i]*surv1[,"survivalD"]+(1-p1[i])
    surv1[,"hazard"]=surv1[,"density"]/surv1[,"survival"]

    surv[[i]]=surv1
  }
  names(surv)=fit$groups.obs$labels


  ### "PD", "PC"
  PD=round(as.vector(p1),3)
  PC=1-PD

  return(list(surv=surv,PD=PD,PC=PC))
}
