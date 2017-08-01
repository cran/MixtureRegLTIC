##########################
### function: split.vars()
##########################
split.vars=function(formula){

  reg=unlist(strsplit(deparse(formula,width.cutoff=500),"~"))[-1]
  reg=unlist(strsplit(reg," ",fixed=TRUE))
  index=which(reg==""|reg=="-"|reg=="-1"|reg=="1"|reg=="+1")
  if(length(index)>0) reg=reg[-index]
  index=which(!is.na(match(reg,"+")))
  if(length(index)>0) reg=reg[-index]

  ### "main", "i2", "i3"
  if(length(reg)==0){
    main=NULL; i2=NULL; i3=NULL
  }else{
    index0=1:length(reg)
    index1=NULL; index2=NULL; index3=NULL
    index=which(!is.na(match(reg,"*")))
    diff.index=c(0,diff(index))
    if(length(index)>0){
      for(i in 1:length(diff.index)){
        if(diff.index[i]==2){
          index3=c(index3,c(index[i]-3,index[i]-1,index[i]+1))
        }else{
          if(i==length(diff.index)){
            index2=c(index2,c(index[i]-1,index[i]+1))
          }else{
            if(diff.index[i+1]!=2) index2=c(index2,c(index[i]-1,index[i]+1))
          }
        }
      }
    }
    index=sort(c(index,index2,index3))
    if(length(index)>0) index1=index0[-index] else index1=index0
    #
    if(length(index1)>0) main=reg[index1] else main=NULL
    if(length(index2)>0) i2=reg[index2] else i2=NULL
    if(length(index3)>0) i3=reg[index3] else i3=NULL

    if(length(main)==0) main=NULL
  }

  ### "y"
  vars=unique(c(main,i2,i3))
  vars.all=all.vars(formula)
  index=match(vars,vars.all)
  if(length(index)>0){
    y=vars.all[-index] 
  }else{
    y=vars.all
  }
  if(length(y)==0) y=NULL
  return(list(y=y,main=main,i2=i2,i3=i3))
}

############################
### function: parameterize()
############################
parameterize=function(formula,eventprobreg,locationreg,scalereg,shapereg,data,shape,var.entry){

  dataset=data

  if(is.null(eventprobreg)) is.reg=list(beta=F) else is.reg=list(beta=T)
  if(is.null(locationreg)) locationreg=~1
  if(is.null(scalereg)) scalereg=~1
  if(is.null(shapereg)) shapereg=~1
  if(is.null(shape)) fix.par=list() else fix.par=list(q=shape)

  vars=split.vars(formula)
  response=vars$y

  #################################
  #  left censored subjects: 2 => 4
  # right censored subjects: 0 => 2
  #################################

  ### "dataset", "var.centime"
  if(length(response)==2){
    time2.add=data.frame(matrix(NA,nrow(dataset),1))
    names(time2.add)=paste(response[1],"2.add",sep="")
    dataset=cbind(dataset,time2.add)
    var.centime=c(response[1],names(time2.add),response[2])
    index=which(dataset[,var.centime[3]]==0) ### 0 ==> 2 (right censored subjects)
    if(length(index)>0){
      dataset[index,var.centime[2]]=NA
      dataset[index,var.centime[3]]=2
    }
    index=which(dataset[,var.centime[3]]==1)
    if(length(index)>0){
      dataset[index,var.centime[2]]=dataset[index,var.centime[1]]
    }
  }
  if(length(response)==3){
    var.centime=response
    index=which(dataset[,var.centime[3]]==2) ### 2 ==> 4 (left censored subjects)
    if(length(index)>0){
      dataset[index,var.centime[1]]=NA
      dataset[index,var.centime[3]]=4
    }
    index=which(dataset[,var.centime[3]]==0) ### 0 ==> 2 (right censored subjects)
    if(length(index)>0){
      dataset[index,var.centime[2]]=NA
      dataset[index,var.centime[3]]=2
    }
  }

  ### "dataset", "var.truntime"
  if(!is.null(var.entry)){
    dataset=cbind(dataset,time2.entry=NA,status.entry=2)
    var.truntime=c(var.entry,"time2.entry","status.entry")
  }else{
    var.truntime=NULL
  }

  ### "is.truncation"
  if(length(var.truntime)==3) is.truncation=T else is.truncation=F

  ### "is.interval"
  if(length(response)==3) is.interval=T else is.interval=F

  ### "covariates"
  covariates=list(names=c(),beta=c(),gamma=c(),alpha=c(),q=c(),
             ibeta=c(),igamma=c(),ialpha=c(),iq=c(),
             i3beta=c(),i3gamma=c(),i3alpha=c(),i3q=c())
  vars.beta=NULL; if(!is.null(eventprobreg)) vars.beta=split.vars(eventprobreg)
  vars.gamma=NULL; if(!is.null(locationreg)) vars.gamma=split.vars(locationreg)
  vars.alpha=NULL; if(!is.null(scalereg)) vars.alpha=split.vars(scalereg)
  vars.q=NULL; if(!is.null(shapereg)) vars.q=split.vars(shapereg)
  vars=unique(c(as.vector(unlist(vars.beta)),as.vector(unlist(vars.gamma)),
              as.vector(unlist(vars.alpha)),as.vector(unlist(vars.q))))
  main.covariates=vars
  if(!is.null(vars)) covariates$names=vars
  if(!is.null(eventprobreg)){
    if(length(vars.beta$main)>0) covariates$beta=match(vars.beta$main,covariates$names)
    if(length(vars.beta$i2)>0) covariates$ibeta=match(vars.beta$i2,covariates$names)
    if(length(vars.beta$i3)>0) covariates$i3beta=match(vars.beta$i3,covariates$names)
  }
  if(!is.null(locationreg)){ 
    if(length(vars.gamma$main)>0) covariates$gamma=match(vars.gamma$main,covariates$names)
    if(length(vars.gamma$i2)>0) covariates$igamma=match(vars.gamma$i2,covariates$names)
    if(length(vars.gamma$i3)>0) covariates$i3gamma=match(vars.gamma$i3,covariates$names)
  }
  if(!is.null(scalereg)){ 
    if(length(vars.alpha$main)>0) covariates$alpha=match(vars.alpha$main,covariates$names)
    if(length(vars.alpha$i2)>0) covariates$ialpha=match(vars.alpha$i2,covariates$names)
    if(length(vars.alpha$i3)>0) covariates$i3alpha=match(vars.alpha$i3,covariates$names)
  }
  if(!is.null(shapereg)){ 
    if(length(vars.q$main)>0) covariates$q=match(vars.q$main,covariates$names)
    if(length(vars.q$i2)>0) covariates$iq=match(vars.q$i2,covariates$names)
    if(length(vars.q$i3)>0) covariates$i3q=match(vars.q$i3,covariates$names)
  }

  #################################################################
  ##### Add two factors and three factors interaction term 
  ##### to dataset and format transformation of variable covariates
  #################################################################
  ### "dataset", "covariates"
  intTerms.Covariates=intCovariates(dataset,covariates)
  dataset=intTerms.Covariates$dataset
  covariates=intTerms.Covariates$covariates
  covariates$main=main.covariates

  num.obs=nrow(dataset)

  ### "x", "covariates"
  x=list();par=list()
  X=data.frame(Intercept=rep(1,num.obs))
  if (!is.null(covariates$names)){
    X=cbind(X,dataset[covariates$names])
  }
  covariates$names=c("Intercept",covariates$names)
  #
  if (is.reg$beta){
    if (is.null(covariates$beta)){
      covariates$beta=1
    }else{
      covariates$beta=c(1,covariates$beta+1)
    }
    x$beta=X[covariates$beta]
    par$beta=rep(0,length(covariates$beta))
  }else{
    covariates$beta=NULL
  }
  #
  if (is.null(covariates$gamma)){
    covariates$gamma=1
  }else{
    covariates$gamma=c(1,covariates$gamma+1)
  }
  x$gamma=X[covariates$gamma]
  par$gamma=rep(0,length(covariates$gamma))
  #
  if (is.null(covariates$alpha)){
    covariates$alpha=1
  }else{
    covariates$alpha=c(1,covariates$alpha+1)
  }
  x$alpha=X[covariates$alpha]
  par$alpha=rep(0,length(covariates$alpha))
  #
  if (is.null(covariates$q)){
    covariates$q=1
  }else{
    covariates$q=c(1,covariates$q+1)
  }
  x$q=X[covariates$q]
  par$q=rep(0,length(covariates$q))

  ### "is.reg"
  is.reg$beta=is.reg$beta
  is.reg$gamma=T
  is.reg$alpha=T
  is.reg$q=T

  ### "minmax.time": Minimum/maximum of survival time
  var=c(var.centime[1:2])
  if(is.truncation){
    if(max(dataset[,var.truntime[1:2]],na.rm=T)>0){
      var=c(var,var.truntime[1:2])
    }
  }

  minmax.time=as.vector(summary(as.vector(as.matrix(dataset[,var])))[c(1,6)])

  return(list(dataset=dataset,covariates=covariates,var.centime=var.centime,var.truntime=var.truntime,is.reg=is.reg,x=x,
         fix.par=fix.par,par=par,is.interval=is.interval,is.truncation=is.truncation,minmax.time=minmax.time))
}


#######################
### function: obsKeep()
#######################
obsKeep=function(dataset,covariates,var.centime,var.truntime){

  num.obs0=dim(dataset)[1]
  use.covariate=NULL
  use.covariate=c(use.covariate,covariates$beta,covariates$ibeta,covariates$i3beta)
  use.covariate=c(use.covariate,covariates$gamma,covariates$igamma,covariates$i3gamma)
  use.covariate=c(use.covariate,covariates$alpha,covariates$ialpha,covariates$i3alpha)
  use.covariate=c(use.covariate,covariates$q,covariates$iq,covariates$i3q)
  use.covariate=unique(use.covariate)

  if(!is.null(use.covariate)) use.covariate=sort(use.covariate)
  if(length(var.truntime)==3){
    var=c(var.centime[3],var.truntime[3],covariates$names[use.covariate])
  }else{
    var=c(var.centime[3],covariates$names[use.covariate])
  }
  
  index=match("Intercept",var);if(length(index)>0) var=var[-index]

  for(k in 1:length(var)){
    index=which(!is.na(dataset[,var[k]]))
    if(length(index)>0) dataset=dataset[index,]
  }

  num.obs=dim(dataset)[1]
  if(num.obs != num.obs0){
    num.diff=num.obs0-num.obs
    print(paste("There is/are ",num.diff," deleted missing subject(s).(",num.obs,"/",num.obs0,")",sep=""))
  }

  return(dataset)
}

#############################
### function: intCovariates()
#############################
####################################################################################################################
##### The dataset of two factors and three factors interaction term and format transformation of variable covariates
####################################################################################################################
intCovariates=function(dataset,covariates){
  ###############
  ### two factors
  ###############
  count=c(covariates$ibeta,covariates$igamma,covariates$ialpha,covariates$iq)
  if(length(count)>=2){
    buffer=NULL;buffer.count=NULL
    for(i in 1:(length(count)/2)){
      buffer1=paste(count[2*i-1],count[2*i],sep="*")
      if(is.null(buffer)){
        buffer=c(buffer,buffer1);buffer.count=c(buffer.count,count[2*i-1],count[2*i])
      }else if(is.na(match(buffer1,buffer))){
        buffer=c(buffer,paste(count[2*i-1],count[2*i],sep="*"))
        buffer.count=c(buffer.count,count[2*i-1],count[2*i])
      }
    }
    count=buffer.count
  }
  idataset=NULL
  if(!is.null(count)){
    idataset=as.data.frame(matrix(NA,dim(dataset)[1],length(count)/2))
    names.icovariates=NULL;labels.icovariates=NULL
    k=1
    for(i in 1:(length(count)/2)){
      names.icovariates=c(names.icovariates,paste(covariates$names[count[2*i-1]],covariates$names[count[2*i]],sep="*"))
      idataset[,k]=dataset[,covariates$names[count[2*i-1]]]*dataset[,covariates$names[count[2*i]]]
      k=k+1
    }
    names(idataset)=names.icovariates
    dataset=cbind(dataset,idataset)

    ### "covariates"
    covariates$names=c(covariates$names,names.icovariates)
    if(!is.null(covariates$ibeta)){
      buffer=NULL
      for(i in 1:(length(covariates$ibeta)/2)){
        location=c(covariates$ibeta[2*i-1],covariates$ibeta[2*i])
        buffer=c(buffer,paste(covariates$names[location[1]],covariates$names[location[2]],sep="*"))
      }
      covariates$ibeta=match(buffer,covariates$names)
      covariates$beta=c(covariates$beta,covariates$ibeta)
    }
    if(!is.null(covariates$igamma)){
      buffer=NULL
      for(i in 1:(length(covariates$igamma)/2)){
        location=c(covariates$igamma[2*i-1],covariates$igamma[2*i])
        buffer=c(buffer,paste(covariates$names[location[1]],covariates$names[location[2]],sep="*"))
      }
      covariates$igamma=match(buffer,covariates$names)
      covariates$gamma=c(covariates$gamma,covariates$igamma)
    }
    if(!is.null(covariates$ialpha)){
      buffer=NULL
      for(i in 1:(length(covariates$ialpha)/2)){
        location=c(covariates$ialpha[2*i-1],covariates$ialpha[2*i])
        buffer=c(buffer,paste(covariates$names[location[1]],covariates$names[location[2]],sep="*"))
      }
      covariates$ialpha=match(buffer,covariates$names)
      covariates$alpha=c(covariates$alpha,covariates$ialpha)
    }
    if(!is.null(covariates$iq)){
      buffer=NULL
      for(i in 1:(length(covariates$iq)/2)){
        location=c(covariates$iq[2*i-1],covariates$iq[2*i])
        buffer=c(buffer,paste(covariates$names[location[1]],covariates$names[location[2]],sep="*"))
      }
      covariates$iq=match(buffer,covariates$names)
      covariates$q=c(covariates$q,covariates$iq)
    }
  }

  #################
  ### three factors
  #################
  count3=c(covariates$i3beta,covariates$i3gamma,covariates$i3alpha,covariates$i3q)
  if(length(count3)>0){
    buffer=NULL;buffer.count3=NULL
    for(i in 1:(length(count3)/3)){
      buffer1=paste(count3[3*i-2],count3[3*i-1],count3[3*i],sep="*")
      if(is.null(buffer)){
        buffer=c(buffer,buffer1);buffer.count3=c(buffer.count3,count3[3*i-2],count3[3*i-1],count3[3*i])
      }else if(is.na(match(buffer1,buffer))){
        buffer=c(buffer,paste(count3[3*i-2],count3[3*i-1],count3[3*i],sep="*"))
        buffer.count3=c(buffer.count3,count3[3*i-2],count3[3*i-1],count3[3*i])
      }
    }
    count3=buffer.count3
  }
  i3dataset=NULL
  if(!is.null(count3)){
    i3dataset=as.data.frame(matrix(NA,dim(dataset)[1],length(count3)/3))
    names.i3covariates=NULL;labels.i3covariates=NULL
    k=1
    for(i in 1:(length(count3)/3)){
      names.i3covariates=c(names.i3covariates,paste(covariates$names[count3[3*i-2]],covariates$names[count3[3*i-1]],covariates$names[count3[3*i]],sep="*"))
      i3dataset[,k]=dataset[,covariates$names[count3[3*i-2]]]*dataset[,covariates$names[count3[3*i-1]]]*dataset[,covariates$names[count3[3*i]]]
      k=k+1
    }
    names(i3dataset)=names.i3covariates
    dataset=cbind(dataset,i3dataset)

    ### "covariates"
    covariates$names=c(covariates$names,names.i3covariates)
    if(!is.null(covariates$i3beta)){
      buffer=NULL
      for(i in 1:(length(covariates$i3beta)/3)){
        location=c(covariates$i3beta[3*i-2],covariates$i3beta[3*i-1],covariates$i3beta[3*i])
        buffer=c(buffer,paste(covariates$names[location[1]],covariates$names[location[2]],covariates$names[location[3]],sep="*"))
      }
      covariates$i3beta=match(buffer,covariates$names)
      covariates$beta=c(covariates$beta,covariates$i3beta)
    }
    if(!is.null(covariates$i3gamma)){
      buffer=NULL
      for(i in 1:(length(covariates$i3gamma)/3)){
        location=c(covariates$i3gamma[3*i-2],covariates$i3gamma[3*i-1],covariates$i3gamma[3*i])
        buffer=c(buffer,paste(covariates$names[location[1]],covariates$names[location[2]],covariates$names[location[3]],sep="*"))
      }
      covariates$i3gamma=match(buffer,covariates$names)
      covariates$gamma=c(covariates$gamma,covariates$i3gamma)
    }
    if(!is.null(covariates$i3alpha)){
      buffer=NULL
      for(i in 1:(length(covariates$i3alpha)/3)){
        location=c(covariates$i3alpha[3*i-2],covariates$i3alpha[3*i-1],covariates$i3alpha[3*i])
        buffer=c(buffer,paste(covariates$names[location[1]],covariates$names[location[2]],covariates$names[location[3]],sep="*"))
      }
      covariates$i3alpha=match(buffer,covariates$names)
      covariates$alpha=c(covariates$alpha,covariates$i3alpha)
    }
    if(!is.null(covariates$i3q)){
      buffer=NULL
      for(i in 1:(length(covariates$i3q)/3)){
        location=c(covariates$i3q[3*i-2],covariates$i3q[3*i-1],covariates$i3q[3*i])
        buffer=c(buffer,paste(covariates$names[location[1]],covariates$names[location[2]],covariates$names[location[3]],sep="*"))
      }
      covariates$i3q=match(buffer,covariates$names)
      covariates$q=c(covariates$q,covariates$i3q)
    }
  }

  covariates=covariates[c("names","beta","gamma","alpha","q")]

  return(list(dataset=dataset,covariates=covariates))
}



