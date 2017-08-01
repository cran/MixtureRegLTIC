#############################
### function: plotResidual()
#############################
plotResidual=function(fit,xlab=NULL,ylab=NULL,main=NULL,col=NULL,lty=NULL,lwd=1,axes=T){

  dataset=fit$res.dataset
  vars=names(dataset)
  main.x=1
  if(length(fit$covariates$main)>0){
    main.x=fit$covariates$main[1]
    if(length(fit$covariates$main)>1){
      for(i in 2:length(fit$covariates$main)){
        main.x=paste(main.x,fit$covariates$main[i],sep="+")
      }
    }
  }
  formu=as.formula(paste("Surv(",vars[1],",",vars[2],",",vars[3],")~",main.x,sep=""))
  var.weight=NULL
  if(fit$is.weight){
    var.weight=vars[length(vars)]
  }
  var.entry=NULL
  if(fit$is.truncation){
    var.entry=vars[4]
    if(all(is.na(dataset[,var.entry]))) var.entry=NULL
  }

  est=NPMLEsurv(formula=formu,data=dataset,var.entry=var.entry,time.origin=0,var.weight=var.weight)

  arg=list(xlab=xlab,ylab=ylab,main=main,col=col,lty=lty,lwd=lwd,axes=axes)

  n.group=length(est$groups.obs$labels)

  if(is.null(arg$xlab)) arg$xlab="Residual"
  if(is.null(arg$ylab)) arg$ylab="Empirical Distribution Function"
  if(is.null(arg$main)){
    arg$main="Empirical distribution curves of residuals"
  }
  if(length(arg$col)!=n.group){
    col=c("black","red","green","blue","cyan","pink","yellow","gray","orange","purple")
    if (n.group > length(col)){
       col=c(col,colors()[-match(col,colors())])
    }
    arg$col=col[1:n.group]
  }
  if(length(arg$lty)!=n.group){
    arg$lty=rep(1,n.group)
  }
  if(is.null(arg$axes)){
    arg$axes=T
  }

  min.x=est$minmax.time[1]
  max.x=est$minmax.time[2]
  x=c(min.x,max.x)
  y=c(0,1)
  plot(x=x,y=y,xlab=arg$xlab,ylab=arg$ylab,type="n",main=arg$main,axes=arg$axes)
  for(i in 1:n.group){
    x=est$surv[[i]][,"time"]
    y=1-est$surv[[i]][,"survivalD"]
    points(x,y,type="s",col=arg$col[i],lty=arg$lty[i],lwd=arg$lwd)
  }

  w=seq(min.x,max.x,0.001)
  buffer=apply(as.matrix(fit$groups.obs$labels),1,strsplit,split=",")
  matrix.cov=matrix(1,length(buffer),1+length(buffer[[1]][[1]]))
  for(k in 1:length(buffer)) matrix.cov[k,]=c(1,as.numeric(buffer[[k]][[1]]))
  if(is.null(fit$fix.par$q)){
    xq=as.matrix(matrix.cov[,fit$covariates$q])%*%as.vector(fit$par$q)
  }else{
    xq=matrix(fit$fix.par$q,n.group)
  }
  q.w=sort(unique(as.vector(xq)))
  q.lty=rep(1:6,each=1,len=length(q.w))
  for(i in 1:length(q.w)){
    if(q.w=="logistic"){
      dist.q=slog(w)
    }else{
      dist.q=sglgd(w,q.w[i])
    }
    points(w,1-dist.q,type="l",col="gray",lty=q.lty[i],lwd=arg$lwd)
  }

  if(length(arg$legend)!=n.group){
    arg$legend=est$groups.obs$labels
  }
  arg$legend=paste(arg$legend," (",est$groups.obs$case.obs," / ",est$groups.obs$total.obs,")",sep="")

  return(arg)
}

############################
### function: residualData()
############################
residualData=function(par,survtime,x){

  ##### "y" and "nobs" (Survival Time and Number of Subjects)
  y=survtime[1:3]
  nobs=dim(y)[1]
  ltau=matrix(0,nobs,1)
  if(length(survtime)==6) ltau=survtime[4]

  ##### "alpha" and "q" (scale and shape)
  alpha=A.fun(x$alpha,par$alpha)
  q=Q.fun(x$q,par$q)

  ##### value of res1.ga, res2.ga and res.tau.ga
  res1.ga=rep(NA,nobs); res2.ga=rep(NA,nobs); res.tau.ga=rep(NA,nobs)
  index=which(y[,3]==1|y[,3]==3|y[,3]==2) # Exact, Interval Censored or Right Censored Data
  if (length(index)>0){
    res1.ga[index]=w.fun(y[index,1],x$gamma[index,],par$gamma,x$alpha[index,],par$alpha)
  }
  index=which(y[,3]==3|y[,3]==4) # Interval or Left Censored Data
  if (length(index)>0){
    res2.ga[index]=w.fun(y[index,2],x$gamma[index,],par$gamma,x$alpha[index,],par$alpha)
  }
  index=which(y[,3]==1) # Exact Data
  if (length(index)>0){
    res2.ga[index]=res1.ga[index]
  }
  if(length(survtime)==6){
    res.tau.ga=w.fun(ltau[,1],x$gamma,par$gamma,x$alpha,par$alpha)
  }
  res.ga.survtime=survtime
  res.ga.survtime[,1]=res1.ga;res.ga.survtime[,2]=res2.ga
  if(length(survtime)==6){
    res.ga.survtime[,4]=res.tau.ga
  }

  return(res.ga.survtime=res.ga.survtime)
}
### example: residualData()
if(F){
  residualData(par,survtime,x)
}
