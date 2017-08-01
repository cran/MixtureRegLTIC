##############################################################
# Find the jump point intervals with Turnbull (1976) algorithm
# corrected by Frydram(1994) and Alioum(1996)
##############################################################
jump.turnbullR=function(survtime,is.truncation,minMax.time){
  ### construct variable 'LRS'
  ### time1 time2 CStatus time3 time4 TStatus
  ###   time index status
  #    time1     3      2 (right censor)
  #    time2     2      4 (left censor)
  #    time1     3      3 (interval censor)
  #    time2     2      3 (interval censor)
  #    time1     1      1 (exact)
  #    time2     2      1 (exact)
  #    time3     2      2 (left truncation)
  #    time4     3      4 (right truncation)
  #    time3     2      3 (interval truncation)
  #    time4     3      3 (interval truncation)
  yc=survtime[,1:3]
  LRS=cbind(time=yc[,1],index=rep(3,length(yc[,1])),status=yc[,3])
  LRS=rbind(LRS,cbind(time=yc[,2],index=rep(2,length(yc[,1])),status=yc[,3]))
  index=which(LRS[,3]==1 & LRS[,2]==3);if(length(index)>0) LRS[index,2]=1
  if(is.truncation){
    yt=survtime[,4:6]
    LRS=rbind(LRS,cbind(time=yt[,1],index=rep(2,length(yt[,1])),status=yt[,3]))
    LRS=rbind(LRS,cbind(time=yt[,2],index=rep(3,length(yt[,1])),status=yt[,3]))
  }
  index=which(!is.na(LRS[,1]));if(length(index)>0) LRS=LRS[index,1:2]
  LRS=rbind(LRS,c(Inf,2));
  if(!is.truncation){
    if(minMax.time[1]<0) LRS=rbind(LRS,c(-Inf,3)) else LRS=rbind(LRS,c(0,3))
  }

  ### Sort the above matrix with second and first columns by the increasing order
  LRS=sort.matrix(M=LRS,col.sort=2,is.increase=T)
  LRS=sort.matrix(M=LRS,col.sort=1,is.increase=T)

  ### Exact data adjustment
  index=which(LRS[,2]==1);if(length(index)>0) LRS[index,2]=3
  buffer=cbind(LRS[,1],c(LRS[-1,1],NA),LRS[,2],c(LRS[-1,2],NA))
  index=which(buffer[,3]>buffer[,4])
  if(is.vector(buffer[index,1:2])){
    jump=as.data.frame(matrix(buffer[index,1:2],1,2))
  }else{
    jump=as.data.frame(buffer[index,1:2])
  }
  names(jump)=c("q","p")
  return(jump)
}
sort.matrix=function(M,col.sort,is.increase=T){
  index=order(M[,col.sort])
  M=M[index,]
  return(M)
}

###########################
### function: turnbullF90()
###########################
turnbullF90=function(survtime,is.truncation=0,weight=NULL,tolerance=1e-6){

  ### minmax.time
  index=which(survtime[,3]==2)
  if(length(index)>0) survtime[index,2]=NA
  index=which(survtime[,3]==4)
  if(length(index)>0) survtime[index,1]=NA
  if (is.truncation==1){
    index=which(survtime[,6]==2)
    if(length(index)>0) survtime[index,5]=NA
    index=which(survtime[,6]==4)
    if(length(index)>0) survtime[index,4]=NA
    if(max(survtime[,c(4,5)],na.rm=T)>0){
      minmax.time=summary(as.vector(as.matrix(survtime[,c(1,2,4,5)])))[c("Min.","Max.")]
    }else{
      minmax.time=summary(as.vector(as.matrix(survtime[,c(1,2)])))[c("Min.","Max.")]
    }
  }else{
    minmax.time=summary(as.vector(as.matrix(survtime[,c(1,2)])))[c("Min.","Max.")]
  }
  minmax.time=as.numeric(minmax.time)

  jump=jump.turnbullR(survtime=survtime,is.truncation=is.truncation,minMax.time=minmax.time)
  ndist=nrow(jump)
  mdist=4

  ### "survtime"
  index=which(is.na(survtime))
  if(length(index)>0) survtime[is.na(survtime)]=999999
  survtime=as.matrix(survtime)
  storage.mode(survtime)="double"

  ### "nsurvtime", "msurvtime", "weight"
  nsurvtime=nrow(survtime)
  msurvtime=ncol(survtime)
  if(is.null(weight)) weight=rep(1,nsurvtime)

  ### "ey", "vary", "skewy", "p0", "p2"
  ey=0; vary=0; skewy=0; p0=0; p2=0

  ### "dist"
  dist=matrix(0,ndist,mdist)
  colnames(dist)=c("q","p","time","survival")
  storage.mode(dist)="double"

  ### Turnbull-Frydman method
  est=.Fortran("turnbulls",ey=as.double(ey),vary=as.double(vary),skewy=as.double(skewy),
               p0=as.double(p0),p2=as.double(p2),survtime=survtime,
               nsurvtime=as.integer(nsurvtime),msurvtime=as.integer(msurvtime),
               istruncation0=as.integer(is.truncation),weight0=as.double(weight),
               dist=dist,ndist=as.integer(ndist),mdist=as.integer(mdist),tolerance0=as.double(tolerance))

  ey=est[1]; vary=est[2]; skewy=est[3]; p0=est[4]; p2=est[5]

  ### "dist"
  dist=as.data.frame(round(est$dist,4))
  if(dist[nrow(dist),"p"]==999999){
    dist[nrow(dist),"p"]=Inf
    dist[nrow(dist),"time"]=round(minmax.time[2],4)
    dist[nrow(dist),"survival"]=dist[nrow(dist)-1,"survival"]
  }
  dist=rbind(data.frame(q=NA,p=NA,time=minmax.time[1],survival=1),dist)
  dist=cbind(dist,survivalD=NA)
  pCure=dist$survival[nrow(dist)]
  dist$survivalD=(dist$survival-pCure)/(1-pCure)
  return(list(dist=dist,E=ey[[1]],Var=vary[[1]],Skew=skewy[[1]]))
}
if(F){
  turnbullF90(survtime,is.truncation=0,weight=NULL,tolerance=1e-6)
}

###############################
### function: turnbullFrydman()
###############################
turnbullFrydman=function(data,var.centime,var.truntime=NULL,covariates=list(names=c(NULL),levels=list(),labels=c(NULL)),
  time.origin=0,var.weight=NULL,tolerance=1e-6){

  time1=Sys.time()

  ### "is.truncation"
  is.truncation=0
  if(length(var.truntime)==3) is.truncation=1

  ### delete missing observations
  var=var.centime[3]
  if(is.truncation) var=c(var,var.truntime[3])
  var=c(var,covariates$names)
  if(!is.null(var.weight)) var=c(var,var.weight)
  index=as.vector(apply(!is.na(data[,var]),1,all))
  if(length(index)>0) data=data[index,]

  ### "minmax.time"
  index=which(data[,var.centime[3]]==2)
  if(length(index)>0) data[index,var.centime[2]]=NA
  index=which(data[,var.centime[3]]==4)
  if(length(index)>0) data[index,var.centime[1]]=NA
  if(length(var.truntime)==3){
    index=which(data[,var.truntime[3]]==2)
    if(length(index)>0) data[index,var.truntime[2]]=NA
    index=which(data[,var.truntime[3]]==4)
    if(length(index)>0) data[index,var.truntime[1]]=NA
  }
  var.survtime=var.centime[1:2]
  if(length(var.truntime)==3){
    if(max(data[,var.truntime[1:2]],na.rm=T)>0){
      var.survtime=c(var.survtime,var.truntime[1:2])
    }
  }
  minmax.time=summary(as.vector(as.matrix(data[,var.survtime])))[c("Min.","Max.")]
  minmax.time=as.vector(minmax.time)

  ### adjust survival time from time origin
  data[,var.centime[1:2]]=data[,var.centime[1:2]]-time.origin
  if(!is.null(var.truntime)){
    data[,var.truntime[1:2]]=data[,var.truntime[1:2]]-time.origin
  }

  ### "groups.obs" (information of groups stratified by covariates)
  groups.obs=list(labels="0",levels=0,levelsTrue=T,obs=rep(0,nrow(data)))
  if(!is.null(covariates$names)){
    groups.obs=groupObs(data=data[,covariates$names],labels=covariates$levels)
  }
  #########################################################
  if(length(covariates$labels)==length(groups.obs$labels)){
    groups.obs$labels=covariates$labels
  }
  #########################################################
  groups.obs$total.obs=rep(0,length(groups.obs$levels))  ### "groups.obs$total.obs"
  groups.obs$case.obs=rep(0,length(groups.obs$levels))   ### "groups.obs$case.obs"
  buffer.data=cbind(data[,var.centime[3]],groups.obs$obs)
  total.obs=table(buffer.data[,2])
  groups.obs$total.obs[match(as.numeric(names(total.obs)),groups.obs$levels)]=total.obs
  index=which(buffer.data[,1] != 2)
  if(length(index)>0) buffer.data=buffer.data[index,]
  case.obs=table(buffer.data[,2])
  groups.obs$case.obs[match(as.numeric(names(case.obs)),groups.obs$levels)]=case.obs

  ### Turnbull-Frydman method
  var.survtime=var.centime
  if(is.truncation) var.survtime=c(var.survtime,var.truntime)
  survtime=data[,var.survtime]
  row.names(survtime)=1:nrow(survtime)
  num.group=length(groups.obs$labels)

  surv=list()
  for (i in 1:num.group){
    index=which(groups.obs$obs==groups.obs$levels[i])
    if(!is.null(var.weight)){
      est=turnbullF90(survtime=survtime[index,],is.truncation=is.truncation,weight=data[index,var.weight],tolerance=tolerance)
    }else{
      est=turnbullF90(survtime=survtime[index,],is.truncation=is.truncation,weight=NULL,tolerance=tolerance)
    }
    est$dist$time=time.origin+est$dist$time
    surv[[i]]=est$dist
  }
  names(surv)=groups.obs$labels

  time2=Sys.time()

  run.time=time2-time1

  run.time=paste(as.numeric(run.time),units(run.time),sep=" ")

  return(list(surv=surv,covariates=covariates,time.origin=time.origin,
                   groups.obs=groups.obs,minmax.time=minmax.time,run.time=run.time))
}
if(F){
  data=cbind(survtime,X)
  var.centime=c("time1","time2","status")
  var.truntime=c("entry","time2.entry","status.entry")
  covariates=list(names="X")
  time.origin=0
  var.weight=NULL
  tolerance=1e-6
  EST=turnbullFrydman(data=data,var.centime=var.centime,var.truntime=var.truntime,covariates=covariates,
                      time.origin=time.origin,var.weight=var.weight,tolerance=tolerance)
}



### R Version ######################################################################################
######################################################
# Turnbull (1976) estimator corrected by Frydram(1994)
######################################################
turnbullR=function(dataset,var.centime,var.truntime=NULL,var.obsweight=NULL,is.graph=F,tolerance=1e-6){
  time1=Sys.time()
  
  ##
  index=which(dataset[,var.centime[3]]==4 & dataset[,var.centime[2]]==0);if(length(index)>0) dataset[index,var.centime[2]]=0.00001
  #print(dataset[,var.centime])
  ##  
  
  ##### Error message
  if(is.null(dataset))return("Dataset Error");if(is.null(var.centime))return("Time Error.")
  if(any(is.na(match(var.centime,names(dataset)))))return("Time Error (Censoring).")
  if(!is.null(var.truntime)){
    if(any(is.na(match(var.truntime,names(dataset)))))return("Time Error (Truncation).")
  }
  if(!is.null(var.obsweight)){
    if(any(is.na(match(var.obsweight,names(dataset))))) return("Weight Error.")
  }

  ##### Keep useful variable and delete missing observations
  dataset=dataset[,c(var.centime,var.truntime,var.obsweight)]
  var=c(var.centime[3],var.obsweight)
  if(length(var.truntime)==3) var=c(var,var.truntime[3])
  for(k in 1:length(var)){
    index=which(!is.na(dataset[,var[k]]))
    if(length(index)>0) dataset=dataset[index,]
  }

  ##### Parameters
  nobs=nrow(dataset)
  is.truncation=F;if(length(var.truntime)==3) is.truncation=T
  survtime=dataset[,c(var.centime,var.truntime)]
  index=which(survtime[,var.centime[3]]==2);if(length(index)>0) survtime[index,var.centime][[2]]=NA
  index=which(survtime[,var.centime[3]]==4);if(length(index)>0) survtime[index,var.centime][[1]]=NA
  if(is.truncation){
    index=which(survtime[,var.truntime][[3]]==2);if(length(index)>0) survtime[index,var.truntime][[2]]=NA
    index=which(survtime[,var.truntime[3]]==4);if(length(index)>0) survtime[index,var.truntime][[1]]=NA
  }
  weight=rep(1,nobs)
  if(!is.null(var.obsweight)) weight=dataset[,var.obsweight]

  time=survtime[,1:2]
  if(is.truncation) time=survtime[,c(1:2,4:5)]
  minMax.time=as.vector(summary(as.vector(as.matrix(time)))[c(1,6)])

  ##### Jump
  jump=jump.turnbullR(survtime,is.truncation,minMax.time)

  ##### Self-consistent estimator
  selfEst=self.turnbullR(survtime,jump,weight,is.truncation,tolerance)

  ##### Graph of Self-consistent estimator
  if(is.graph){
    Est=selfEst[,c("time","survival")]
    Est=rbind(c(Est[1,"time"],1),Est)  
    if(Est[nrow(Est),"time"]==Inf & nrow(Est)>1){
      Est[nrow(Est),"time"]=minMax.time[2]*(1+0.05)
      Est[nrow(Est),"survival"]=Est[nrow(Est)-1,"survival"]
    }
    plot(Est[,"time"],Est[,"survival"],type="s",ylim=c(0,1),xlab="Time",ylab="Probability",main="Survival Curve Estimator\n (Turnbull Method)")
  }
  time2=Sys.time()
  time.diff=time2-time1
#  print(time.diff)
  return(selfEst)
}
##############################################################
# Find the jump point intervals with Turnbull (1976) algorithm
# corrected by Frydram(1994) and Alioum(1996)
##############################################################
jump.turnbullR=function(survtime,is.truncation,minMax.time){
  ### construct variable 'LRS'
  ### time1 time2 CStatus time3 time4 TStatus
  ###   time index status
  #    time1     3      2 (right censor) 
  #    time2     2      4 (left censor)
  #    time1     3      3 (interval censor)
  #    time2     2      3 (interval censor)
  #    time1     1      1 (exact)
  #    time2     2      1 (exact)
  #    time3     2      2 (left truncation) 
  #    time4     3      4 (right truncation)
  #    time3     2      3 (interval truncation)
  #    time4     3      3 (interval truncation)
  yc=survtime[,1:3]
  LRS=cbind(time=yc[,1],index=rep(3,length(yc[,1])),status=yc[,3])
  LRS=rbind(LRS,cbind(time=yc[,2],index=rep(2,length(yc[,1])),status=yc[,3]))
  index=which(LRS[,3]==1 & LRS[,2]==3);if(length(index)>0) LRS[index,2]=1
  if(is.truncation){
    yt=survtime[,4:6]
    LRS=rbind(LRS,cbind(time=yt[,1],index=rep(2,length(yt[,1])),status=yt[,3]))
    LRS=rbind(LRS,cbind(time=yt[,2],index=rep(3,length(yt[,1])),status=yt[,3]))
  }
  index=which(!is.na(LRS[,1]));if(length(index)>0) LRS=LRS[index,1:2]
  LRS=rbind(LRS,c(Inf,2));
  if(!is.truncation){
    if(minMax.time[1]<0) LRS=rbind(LRS,c(-Inf,3)) else LRS=rbind(LRS,c(0,3))
  }

  ### Sort the above matrix with second and first columns by the increasing order
  LRS=sort.matrix(M=LRS,col.sort=2,is.increase=T)
  LRS=sort.matrix(M=LRS,col.sort=1,is.increase=T)

  ### Exact data adjustment
  index=which(LRS[,2]==1);if(length(index)>0) LRS[index,2]=3
  buffer=cbind(LRS[,1],c(LRS[-1,1],NA),LRS[,2],c(LRS[-1,2],NA))
  index=which(buffer[,3]>buffer[,4])
  if(is.vector(buffer[index,1:2])){
    jump=as.data.frame(matrix(buffer[index,1:2],1,2))
  }else{
    jump=as.data.frame(buffer[index,1:2])
  }
  names(jump)=c("q","p")
  return(jump)
}
sort.matrix=function(M,col.sort,is.increase=T){
  index=order(M[,col.sort])
  M=M[index,]
  return(M)
}
######################################################
# Find the mass of jump intervals with Turnbull (1976)
# self-consistent algorithm
######################################################
self.turnbullR=function(survtime,jump,weight=1,is.truncation=F,tolerance=1e-4){
  alpha=AlphaBeta(survtime[,1:3],jump)
  if(is.truncation) beta=AlphaBeta(survtime[,4:6],jump)

  max.iteration=Inf
  njump=dim(jump)[1]
  mjump=rep(NA,njump)
  mjump0=mjump;mjump0[]=1/njump

  iteration=0;diff=Inf
  while(iteration<max.iteration & diff>tolerance){
    iteration=iteration+1
    muij=t(t(alpha)*mjump0)*weight/as.vector(alpha%*%mjump0)
    Eij=muij
    if(is.truncation){
      nuij=t(t((1-beta))*mjump0)*weight/as.vector(beta%*%mjump0)
      Eij=muij+nuij
    }
    mjump[]=apply(Eij,2,sum)/sum(Eij)
    diff=max(abs(mjump-mjump0))
    mjump0=mjump
  }
  survival=rep(0,njump)
  for(j in 2:njump){
    survival[j-1]=sum(mjump[j:njump])
  }

  value.round=1
  for(k in 1:6){
    if(tolerance*10^k >=1){
      value.round=k
      break
    }
  }
  survival=round(survival,value.round)

  return(cbind(jump,time=jump[,2],survival))
}
AlphaBeta=function(time,jump){
  AB=matrix(0,nrow(time),nrow(jump))
  # observation time
  timel=time[,1];timer=time[,2]
  index=which(is.na(timel));if(length(index)>0) timel[index]=-Inf
  index=which(is.na(timer));if(length(index)>0) timer[index]=Inf
  # jump time
  jumpl=matrix(1,nrow(jump),nrow(time))
  jumpl[]=jump[,1];jumpl=t(jumpl)
  jumpr=matrix(1,nrow(jump),nrow(time))
  jumpr[]=jump[,2];jumpr=t(jumpr)
  # 
  AB[timel<=jumpl & jumpr<=timer]=1
  AB[timel==jumpl & timel==jumpr & timel!=timer]=0 # for exact time
  colnames(AB)=paste(jump[,1],jump[,2],sep="-")
  rownames(AB)=paste(timel,timer,sep="-")
  return(AB)
}
#######################################
# Find the conditional distribution and 
# it's mean, variance and skewness 
#######################################
condition.turnbullR=function(turnbullEst,is.graph=F,restime=NULL){
  ### Distribution of susceptibility
  SDEst=turnbullEst[,c("time","survival")]
  if(!is.null(restime)){
    index=which(SDEst$time>restime)
    if(length(index)) SDEst=SDEst[-index,]
  }
  if(SDEst[1,"time"]==-Inf) SDEst=SDEst[-1,]
  if(SDEst[nrow(SDEst),"time"]==Inf) SDEst=SDEst[-nrow(SDEst),]
  pcure=SDEst[nrow(SDEst),"survival"]
  SDEst[,"survival"]=(SDEst[,"survival"]-pcure)/(1-pcure)

  ### Mean and variance of susceptibility
  zero.time=!(length(which(SDEst[,"time"]==0))==0)
  if(!zero.time){
    SDEst=rbind(SDEst,c(0,NA));SDEst=SDEst[order(SDEst[,"time"]),]
    index=which(SDEst[,"time"]==0)
    if(index==1){ 
      SDEst[index,"survival"]=1
    }else if(index>1){ 
      SDEst[index,"survival"]=SDEst[index-1,"survival"]
    }
  }
  ### EY, EY2 and EY3
  EY=0;EY2=0;EY3=0
  index=which(SDEst[,"time"]>=0)
  if(length(index)>0){
    jumpSur=SDEst[index,c("time","survival")]
    time=jumpSur[,"time"];Sur=jumpSur[-nrow(jumpSur),"survival"]
    EY=EY+sum(diff(time)*Sur)
    EY2=EY2+sum(diff(time^2)*Sur)
    EY3=EY3+sum(diff(time^3)*Sur)
  }
  index=which(SDEst[,"time"]<=0)
  if(length(index)>0){
    jumpSur=SDEst[index,c("time","survival")]
    time=jumpSur[,"time"];Sur=jumpSur[-nrow(jumpSur),"survival"]
    EY=EY-sum(diff(time)*(1-Sur))
    EY2=EY2-sum(diff(time^2)*(1-Sur))
    EY3=EY3-sum(diff(time^3)*(1-Sur))
  }
  ### VarY and SkewY
  VarY=EY2-(EY^2)
  SkewY=EY3-3*EY*VarY-(EY^3)

  if(!zero.time){
    index=which(SDEst[,"time"]==0)
    if(length(index)>0) SDEst=SDEst[-index,]
  }

#  ##### Graph of Self-consistent estimator
#  if(is.graph){
#    Est=SDEst[,c("time","survival")]
#    Est=rbind(c(Est[1,"time"],1),Est)  
#    if(Est[nrow(Est),"time"]==Inf & nrow(Est)>1){
#      Est[nrow(Est),"time"]=minMax.time[2]*(1+0.05)
#      Est[nrow(Est),"survival"]=Est[nrow(Est)-1,"survival"]
#    }
#    plot(Est[,"time"],Est[,"survival"],type="s",ylim=c(0,1),xlab="Time",ylab="Probability",main="Conditional Survival Curve Estimator\n (Turnbull Method)")
#  }

  return(list(SDEst=SDEst,EY=round(EY,10),VarY=round(VarY,10),SkewY=round(SkewY,10)))
}

############
### Example:
############
if(F){
  ### Ex1:
  dataset=read.table("data_turnbull.txt",header=T,sep="\t")
  var.centime=names(dataset)[1:3]
  var.truntime=NULL
  var.obsweight= NULL
  is.graph=F
  tolerance=1e-6
  turnbullEst=turnbullR(dataset,var.centime,var.truntime,var.obsweight,is.graph,tolerance)
  turnbullEst
#  cond.turnbullEst=condition.turnbullR(turnbullEst,is.graph,restime=NULL)
#  cond.turnbullEst
#  cond.turnbullEst=condition.turnbullR(turnbullEst,is.graph,restime=20000)
#  cond.turnbullEst

  ### Ex2:
  dataset=read.table("ICen.dat",header=T,sep="")
  var.centime=names(dataset)[1:3]
  var.truntime=NULL
  var.obsweight=names(dataset)[6]
  is.graph=F
  tolerance=1e-6
  turnbullEst=turnbullR(dataset,var.centime,var.truntime,var.obsweight,is.graph,tolerance)
  turnbullEst
#  cond.turnbullEst=condition.turnbullR(turnbullEst,is.graph,restime=NULL)
#  cond.turnbullEst

  ### Ex3:
  dataset=read.table("LTrunICen.dat",header=T,sep=" ")
  var.centime=names(dataset)[1:3]
  var.truntime=names(dataset)[4:6]
  var.obsweight=names(dataset)[9]
  is.graph=F
  tolerance=1e-6
  turnbullEst=turnbullR(dataset,var.centime,var.truntime,var.obsweight,is.graph,tolerance)
  turnbullEst
#  cond.turnbullEst=condition.turnbullR(turnbullEst,is.graph,restime=NULL)
#  cond.turnbullEst

  ### Ex4:
  dataset=read.table("test.dat",header=T,sep=" ")
  var.centime=names(dataset)[1:3]
  var.truntime=names(dataset)[4:6]
  var.obsweight=NULL
  is.graph=F
  tolerance=1e-6
  turnbullEst=turnbullR(dataset,var.centime,var.truntime,var.obsweight,is.graph,tolerance)
  turnbullEst
#  cond.turnbullEst=condition.turnbullR(turnbullEst,is.graph,restime=NULL)
#  cond.turnbullEst

}

