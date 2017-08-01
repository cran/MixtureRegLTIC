##############################
### function: plotNPMLEsurv()
##############################
plotNPMLEsurv=function(est,dist="overall",curve="survival",type="s",
  xlab=NULL,ylab=NULL,main=NULL,col=NULL,lty=NULL,lwd=1,axes=T){

  arg=list(type=type,xlab=xlab,ylab=ylab,main=main,col=col,lty=lty,lwd=lwd,axes=axes)

  n.group=length(est$groups.obs$labels)

  if(is.null(arg$xlab)){
    arg$xlab="Time"
  }
  if(is.null(arg$ylab)){
    arg$ylab="Survival Probability"
    if(curve=="event") arg$ylab="Event Probability"
  }
  if(is.null(arg$main)){
    arg$main="Survival curves"
    if(curve=="event") arg$main="Event curves"
    if(dist=="cond"){
      arg$main="Conditional survival curves \n in logarithmic scale"
      if(curve=="event") arg$main="Conditional event curves \n in logarithmic scale"
    }
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

  if(dist=="cond"){
    if(est$time.origin>0){
      arg$xlab=paste("log(",arg$xlab," - ",est$time.origin,")",sep="")
    }else{
      arg$xlab=paste("log(",arg$xlab,")",sep="")
    }
    x=log(est$minmax.time-est$time.origin)
    y=c(0,1)
    plot(x=x,y=y,type="n",xlab=arg$xlab,ylab=arg$ylab,main=arg$main,axes=arg$axes)
    for(i in 1:n.group){
      x=log(est$surv[[i]][,"time"]-est$time.origin)
      if(curve=="event"){
        y=1-est$surv[[i]][,"survivalD"]
      }else{
        y=est$surv[[i]][,"survivalD"]
      }
      points(x,y,type=arg$type,col=arg$col[i],lty=arg$lty[i],lwd=arg$lwd)
    }
  }else{
    x=est$minmax.time
    y=c(0,1)
    plot(x=x,y=y,xlab=arg$xlab,ylab=arg$ylab,type="n",main=arg$main,axes=arg$axes)
    for(i in 1:n.group){
      x=est$surv[[i]][,"time"]
      if(curve=="event"){
        y=1-est$surv[[i]][,"survival"]
      }else{
        y=est$surv[[i]][,"survival"]
      }
      points(x,y,type=arg$type,col=arg$col[i],lty=arg$lty[i],lwd=arg$lwd)
    }
  }

  if(length(arg$legend)!=n.group){
    arg$legend=est$groups.obs$labels
  }
  arg$legend=paste(arg$legend," (",est$groups.obs$case.obs," / ",est$groups.obs$total.obs,")",sep="")

  return(arg)
}
