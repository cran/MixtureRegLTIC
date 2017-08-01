#############################
### function: printMixture()
#############################
printMixture=function(fit,digits=3,file=NULL){

  est=fit

  names.var=est$covariates$names
  is.truncation=est$is.truncation
  is.interval=est$is.interval
  n.beta=length(est$covariates$beta)
  n.gamma=length(est$covariates$gamma)
  n.alpha=length(est$covariates$alpha)
  n.q=length(est$covariates$q)

  if(is.null(file)){

    cat("\n")
    cat("Call: "); print(est$call)
    cat("\n")

    if(length(est$par$beta)==0){
      cat("The AFT Location-Scale Mixture Regression Model for ")
    }else{
      cat("The Logistic-AFT Location-Scale Mixture Regression Model for ")
    }
    if(is.truncation){
      if(is.interval){
        cat("Left-Truncated and Interval-Censored Data \n")
      }else{
        cat("Left-Truncated and Right-Censored Data \n")
      }
    }else{
      if(is.interval){
        cat("Interval-Censored Data \n")
      }else{
        cat("Right-Censored Data \n")
      }
    }
    cat("\n")

    empty="    "
    ### "rownames.outTable", "colnames.outTable"
    colnames.outTable=c("EST","STD","95% LCL","95% UCL","P-Value","Exp(EST)","Gradient")
    rownames.outTable=NULL
    if(n.beta>0){
      rownames.outTable=c(rownames.outTable,"Event Probability Sub-model ")
      loc=est$covariates$beta
      for (i in 1:length(loc)){
        rownames.outTable=c(rownames.outTable,paste(empty,names.var[loc[i]],sep=""))
      }
    }
    rownames.outTable=c(rownames.outTable,"AFT Location-Scale Sub-model")
    rownames.outTable=c(rownames.outTable,"  Location Part")
    loc=est$covariates$gamma
    for (i in 1:length(loc)){
      rownames.outTable=c(rownames.outTable,paste(empty,names.var[loc[i]],sep=""))
    }
    rownames.outTable=c(rownames.outTable,"  Scale Part")
    loc=est$covariates$alpha
    for (i in 1:length(loc)){
      rownames.outTable=c(rownames.outTable,paste(empty,names.var[loc[i]],sep=""))
    }
    if(is.null(est$fix.par$q)){
      if(n.q==1){
        rownames.outTable=c(rownames.outTable,"  Shape Parameter")
      }else{
        rownames.outTable=c(rownames.outTable,"  Shape Part")
        loc=est$covariates$q
        for (i in 1:length(loc)){
          rownames.outTable=c(rownames.outTable,paste(empty,names.var[loc[i]],sep=""))
        }
      }
    }else{
      if(est$fix.par$q=="logistic"){
        rownames.outTable=c(rownames.outTable,"  Log-Logistic Distribution")
      }else{
        rownames.outTable=c(rownames.outTable,"  Shape Parameter")
      }
    }

    ### "outTable"
    outTable=as.data.frame(matrix(NA,length(rownames.outTable),length(colnames.outTable)))
    outTable=matrix(NA,length(rownames.outTable),length(colnames.outTable))
#    outTable=as.data.frame(outTable)
    colnames(outTable)=colnames.outTable
    row.names(outTable)=rownames.outTable

    k=0
    ##### Logistic Regression (Event Probability)
    j=0
    if(n.beta>0){
      j=j+1
      loc=est$covariates$beta
      for (i in 1:length(loc)){
        j=j+1
        k=k+1
        parest=round(est$parest[k],digits)
        std=round(est$std[k],digits)
        LCL=round(parest-qnorm(1-0.025)*std,digits)
        UCL=round(parest+qnorm(1-0.025)*std,digits)
        pvalue=round(est$pvalue[k],digits)
        if(i==1) exp.est=NA else exp.est=round(exp(est$parest[k]),digits)
        grad=est$gradient[k]
        outTable[j,]=c(parest,std,LCL,UCL,pvalue,exp.est,grad)
      }
    }

    ### Location Part
    j=j+2
    loc=est$covariates$gamma
    for (i in 1:length(loc)){
      j=j+1
      k=k+1
      parest=round(est$parest[k],digits)
      std=round(est$std[k],digits)
      LCL=round(parest-qnorm(1-0.025)*std,digits)
      UCL=round(parest+qnorm(1-0.025)*std,digits)
      pvalue=round(est$pvalue[k],digits)
      exp.est=NA
      grad=est$gradient[k]
      outTable[j,]=c(parest,std,LCL,UCL,pvalue,exp.est,grad)
    }

    ### Scale Part
    j=j+1
    loc=est$covariates$alpha
    for (i in 1:length(loc)){
      j=j+1
      k=k+1
      parest=round(est$parest[k],digits)
      std=round(est$std[k],digits)
      LCL=round(parest-qnorm(1-0.025)*std,digits)
      UCL=round(parest+qnorm(1-0.025)*std,digits)
      pvalue=round(est$pvalue[k],digits)
      exp.est=NA
      grad=est$gradient[k]
      outTable[j,]=c(parest,std,LCL,UCL,pvalue,exp.est,grad)
    }

    ### Shape Parameter/Part
    if(is.null(est$fix.par$q) & n.q!=1) j=j+1
    loc=est$covariates$q
    for (i in 1:length(loc)){
      j=j+1
      k=k+1
      if(!is.null(est$fix.par$q)){
        if(est$par$q=="logistic"){
          parest=est$par$q
        }else{
          parest=round(est$par$q,digits)
        }
        outTable[j,1]=parest
      }else{
        parest=round(est$parest[k],digits)
        std=round(est$std[k],digits)
        LCL=round(parest-qnorm(1-0.025)*std,digits)
        UCL=round(parest+qnorm(1-0.025)*std,digits)
        pvalue=round(est$pvalue[k],digits)
        exp.est=NA
        grad=est$gradient[k]
        outTable[j,]=c(parest,std,LCL,UCL,pvalue,exp.est,grad)
      }
    }
    print(outTable[,-6],na.print="",quote=F)
    cat("\n")
    cat(paste("Log-Likelihood = ",round(est$LLF,digits),"\n",sep=""))
    cat(paste("AIC = ",round(est$AIC,digits),"\n",sep=""))

    if(is.null(est$COV)){
      cat("warning: Hessian Matrix is Singular !! \n")
    }else if(est$convergence!=0){
      cat(paste("warning: No Convergence !! \n",sep=","))
    }
    cat("\n")

  }else{
    write("",file=file,append=F)
    if(length(est$par$beta)==0){
      write(paste("","","The AFT Location-Scale Mixture Regression Model",sep=","),file=file,append=T)
    }else{
      write(paste("","","The Logistic-AFT Location-Scale Mixture Regression Model",sep=","),file=file,append=T)
    }
    if(is.truncation){
      if(is.interval){
        write(paste("","","            for Left-Truncated and Interval-Censored Data",sep=","),file=file,append=T)
      }else{
        write(paste("","","             for Left-Truncated and Right-Censored Data",sep=","),file=file,append=T)
      }
    }else{
      if(is.interval){
        write(paste("","","                        for Interval-Censored Data",sep=","),file=file,append=T)
      }else{
        write(paste("","","                          for Right-Censored Data",sep=","),file=file,append=T)
      }
    }
    write("",file=file,append=T)
    write(paste("","","","EST","STD","95% LCL","95% UCL","P-Value","Gradient",sep=","),file=file,append=T)

    k=0
    #############################################
    ##### Logistic Regression (Event Probability)
    #############################################
    if(n.beta>0){
      write(paste("","Event Probability Sub-model",sep=","),file=file,append=T)
      loc=est$covariates$beta
      for (i in 1:length(loc)){
        k=k+1
        parest=round(est$parest[k],digits)
        std=round(est$std[k],digits)
        LCL=round(parest-qnorm(1-0.025)*std,digits)
        UCL=round(parest+qnorm(1-0.025)*std,digits)
        pvalue=ifelse(round(est$pvalue[k],digits)==0,"      <0.001",round(est$pvalue[k],digits))
        grad=ifelse(round(est$gradient[k],digits=10)==0,"<1.0E-10",round(est$gradient[k],digits=10))
        if (i==1){
          write(paste("","",names.var[loc[i]],parest,std,LCL,UCL,pvalue,grad,sep=","),file=file,append=T)
        }else{
          write(paste("","",names.var[loc[i]],parest,std,LCL,UCL,pvalue,grad,sep=","),file=file,append=T)
        }
      }
    }

    ##### AFT Regression (Event Time Distribution)
    write(paste("","AFT Location-Scale Sub-model",sep=","),file=file,append=T)

    ### Location Part
    write(paste("","  Location Part",sep=","),file=file,append=T)
    loc=est$covariates$gamma
    for (i in 1:length(loc)){
      k=k+1
      parest=round(est$parest[k],digits)
      std=round(est$std[k],digits)
      LCL=round(parest-qnorm(1-0.025)*std,digits)
      UCL=round(parest+qnorm(1-0.025)*std,digits)
      pvalue=ifelse(round(est$pvalue[k],digits)==0, "      <0.001",round(est$pvalue[k],digits))
      grad=ifelse(round(est$gradient[k],digits=10)==0,"<1.0E-10",round(est$gradient[k],digits=10))
      if (i==1) write(paste("","",names.var[loc[i]],parest,std,LCL,UCL,pvalue,grad,sep=","),file=file,append=T)
      else write(paste("","",names.var[loc[i]],parest,std,LCL,UCL,pvalue,grad,sep=","),file=file,append=T)
    }

    ### Scale Part
    write(paste("","  Scale Part",sep=","),file=file,append=T)
    loc=est$covariates$alpha
    for (i in 1:length(loc)){
      k=k+1
      parest=round(est$parest[k],digits)
      std=round(est$std[k],digits)
      LCL=round(parest-qnorm(1-0.025)*std,digits)
      UCL=round(parest+qnorm(1-0.025)*std,digits)
      pvalue=ifelse(round(est$pvalue[k],digits)==0,"      <0.001",round(est$pvalue[k],digits))
      grad=ifelse(round(est$gradient[k],digits=10)==0,"<1.0E-10",round(est$gradient[k],digits=10))
      if (i==1) write(paste("","",names.var[loc[i]],parest,std,LCL,UCL,pvalue,grad,sep=","),file=file,append=T)
      else write(paste("","",names.var[loc[i]],parest,std,LCL,UCL,pvalue,grad,sep=","),file=file,append=T)
    }

    ### Shape Component
    if(is.null(est$fix.par$q)){
      if(n.q==1){
        write(paste("","  Shape Parameter",sep=","),file=file,append=T)
      }else{
        write(paste("","  Shape Part",sep=","),file=file,append=T)
      }
    }else{
      if(est$fix.par$q=="logistic"){
        write(paste("","  Log-Logistic Distribution",sep=","),file=file,append=T)
      }else{
        write(paste("","Shape Parameter",sep=","),file=file,append=T)
      }
    }
    loc=est$covariates$q
    for (i in 1:length(loc)){
      if(is.null(est$fix.par$q)){
        k=k+1
        parest=round(est$parest[k],digits)
        std=round(est$std[k],digits)
        LCL=round(parest-qnorm(1-0.025)*std,digits)
        UCL=round(parest+qnorm(1-0.025)*std,digits)
        pvalue=ifelse(round(est$pvalue[k],digits)==0,"      <0.001",round(est$pvalue[k],digits))
        grad=ifelse(round(est$gradient[k],digits=10)==0,"<1.0E-10",round(est$gradient[k],digits=10))
        if (i==1 & length(loc)==1){
          write(paste("","","",parest,std,LCL,UCL,pvalue,"",grad,sep=","),file=file,append=T)
        }else if (i==1){
          write(paste("","",names.var[loc[i]],parest,std,LCL,UCL,pvalue,grad,sep=","),file=file,append=T)
        }
        else write(paste("","",names.var[loc[i]],parest,std,LCL,UCL,pvalue,grad,sep=","),file=file,append=T)
      }else{
        if(length(loc)==1){
          if(est$par$q[i]!="logistic"){
            write(paste("","","",round(est$par$q[i],digits),sep=","),file=file,append=T)
          }
        }else{
          write(paste("","",names.var[loc[i]],"",sep=","),file=file,append=T)
        }
      }
    }
    write("",file=file,append=T)
    write(paste(",Log-Likelihood = ",round(est$LLF,digits),sep=""),file=file,append=T)
    write(paste(",AIC = ",round(est$AIC,digits),sep=""),file=file,append=T)

    if(is.null(est$COV)){
      write("warning: Hessian Matrix is Singular !!",file=file,append=T)
    }else if(est$convergence!=0){
      write(paste("","warning: No Convergence !!",sep=","),file=file,append=T)
    }
  }
}
