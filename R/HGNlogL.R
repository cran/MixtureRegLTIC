###################
### function: LLF()
###################
LLF=function(par.nfix,par.all,is.fix.par,survtime,x,is.reg,weight,mtype){
  LLF=HGNlogL(par.nfix,par.all,is.fix.par,survtime,x,is.reg,weight,mtype,outtype="LLF")
  return(LLF)
}

####################
### function: GLLF()
####################
GLLF=function(par.nfix,par.all,is.fix.par,survtime,x,is.reg,weight,mtype){
  GLLF=HGNlogL(par.nfix,par.all,is.fix.par,survtime,x,is.reg,weight,mtype,outtype="GLLF")
  return(GLLF)
}

####################
### function: HLLF()
####################
HLLF=function(par.nfix,par.all,is.fix.par,survtime,x,is.reg,weight,mtype){
  HLLF=HGNlogL(par.nfix,par.all,is.fix.par,survtime,x,is.reg,weight,mtype,outtype="HLLF")
  return(HLLF)
}

#######################
### function: HGNlogL()
#######################
HGNlogL=function(par.nfix,par.all,is.fix.par,survtime,x,is.reg,weight,mtype,outtype="LLF"){
  HGNlogL=0
  if (length(survtime)==6){
    # Calculate truncation contribution
    HGNlogL=HGNlogL-HGNlogLC(par.nfix,par.all,is.fix.par,survtime[4:6],x,is.reg,weight,mtype,outtype)
  }
  # Calculate interval censorship contribution
  HGNlogL=HGNlogL+HGNlogLC(par.nfix,par.all,is.fix.par,survtime[1:3],x,is.reg,weight,mtype,outtype)
  return(HGNlogL)
}

########################
### function: HGNlogLC()
########################
HGNlogLC=function(par.nfix,par.all,is.fix.par,survtime,x,is.reg,weight,mtype,outtype="LLF"){

  ##### "par.all" => "par"
  par=par.all

  ##### "par.nfix" => "par"
  start=1
  if (is.reg$beta){
    len=length(par$beta); par$beta=par.nfix[seq(start,start+len-1)]; start=start+len
  }
  len=length(par$gamma); par$gamma=par.nfix[seq(start,start+len-1)]; start=start+len
  len=length(par$alpha); par$alpha=par.nfix[seq(start,start+len-1)]; start=start+len
  if (!is.fix.par$q){
    len=length(par$q); par$q=par.nfix[seq(start,start+len-1)]; start=start+len
  }

  ##### "y" and "nobs" (Survival Time and Number of Subjects)
  y=survtime;nobs=dim(y)[1]

  ##### Mixture Type (The Type of Distribution Form for Each Subject)
  ip2=rep(1,nobs) # "two components with cure"
  index=which(mtype==1) # "one component"
  if (length(index)>0){ip2[index]=0}
  ip1=ip2 # used for derivative with beta

  #######################################
  ##### "alpha" and "q" (scale and shape)
  #######################################
  alpha=A.fun(x$alpha,par$alpha)
  q=Q.fun(x$q,par$q)

  ##########################################
  ##### Functions of "p1" & "p2=1-p1"
  ##### and It's First and Second Derivative
  ##########################################
  ##### "p1","p2" (function)
  p1=rep(1,nobs)
  if(is.reg$beta){
    p1=p1.fun(x,par,is.reg)
    index=which(mtype==1|mtype==3);if(length(index)>0) p1[index]=1
  }
  p2=rep(0,nobs);p2=1-p1
  
  if(outtype=="GLLF"|outtype=="HLLF"){
    ##### "gp1","gp2" (first derivative)
    gp1=list(beta=rep(0,nobs))
    if (is.reg$beta){ 
      gp1$beta=p1b.fun(x,par,is.reg)
      index=which(mtype==1|mtype==3);if(length(index)>0) gp1$beta[index]=0
    }
    gp2=list(beta=rep(0,nobs));gp2$beta=-gp1$beta

    if(outtype=="HLLF"){
      ##### "hp1", "hp2" (second derivative)
      hp1=list(beta=rep(0,nobs))
      if (is.reg$beta){ 
        hp1$beta=p1bb.fun(x,par,is.reg)
        index=which(mtype==1|mtype==3);if (length(index)>0) hp1$beta[index]=0
      }
      hp2=list(beta=rep(0,nobs));hp2$beta=-hp1$beta
    }
  }

  ################
  ##### Value of w
  ################
  ##### "w1","w2"
  w1=rep(NA,nobs);w2=rep(NA,nobs)
  index=which(y[,3]==1|y[,3]==3|y[,3]==2) # Exact, Interval Censored or Right Censored Data
  if (length(index)>0){w1[index]=w.fun(y[index,1],x$gamma[index,],par$gamma,x$alpha[index,],par$alpha)}
  index=which(y[,3]==3|y[,3]==4) # Interval or Left Censored Data
  if (length(index)>0){w2[index]=w.fun(y[index,2],x$gamma[index,],par$gamma,x$alpha[index,],par$alpha)}

  ##########################################
  ##### Density and Survival Functions of w
  ##### and It's First and Second Derivative
  ##########################################
  ##### Density and Survival Functions of w
  ### "fw1","Sw1","fw2","Sw2"
  fw1=rep(NA,nobs);Sw1=rep(NA,nobs)
  fw2=rep(NA,nobs);Sw2=rep(NA,nobs)
  index=which(y[3]==1|y[3]==3|y[3]==2) # Exact, Interval Censored or Right Censored Data
  if (length(index)>0){
    fw1[index]=dw.fun(w1[index],q[index])
    Sw1[index]=Sw.fun(w1[index],q[index])
  }
  index=which(y[3]==3|y[3]==4) ### Interval Censored or Left Censored Data
  if (length(index)>0){
    fw2[index]=dw.fun(w2[index],q[index])
    Sw2[index]=Sw.fun(w2[index],q[index])
  }

  if(outtype=="GLLF"|outtype=="HLLF"){
    ##### First Derivative of Density and Survival Functions of w
    ### "gfw1","gSw1","gfw2","gSw2"
    gfw1=list(gamma=rep(NA,nobs),alpha=rep(NA,nobs),q=rep(NA,nobs))
    gSw1=list(gamma=rep(NA,nobs),alpha=rep(NA,nobs),q=rep(NA,nobs))
    gfw2=list(gamma=rep(NA,nobs),alpha=rep(NA,nobs),q=rep(NA,nobs))
    gSw2=list(gamma=rep(NA,nobs),alpha=rep(NA,nobs),q=rep(NA,nobs))

    index=which(y[3]==1|y[3]==3|y[3]==2) # Exact, Interval Censored or Right Censored Data
    if (length(index)>0){
      gfw1$gamma[index]=dwg.fun(w1[index],q[index],y[index,1],x$gamma[index,],par$gamma,x$alpha[index,],par$alpha)
      gfw1$alpha[index]=dwa.fun(w1[index],q[index],y[index,1],x$gamma[index,],par$gamma,x$alpha[index,],par$alpha)
      if(!is.fix.par$q){
        gfw1$q[index]=dwq.fun(w1[index],q[index])
      }
      gSw1$gamma[index]=Swg.fun(w1[index],q[index],y[index,1],x$gamma[index,],par$gamma,x$alpha[index,],par$alpha)
      gSw1$alpha[index]=Swa.fun(w1[index],q[index],y[index,1],x$gamma[index,],par$gamma,x$alpha[index,],par$alpha)
      if(!is.fix.par$q){gSw1$q[index]=Swq.fun(w1[index],q[index])}
    }
    index=which(y[3]==3|y[3]==4) # Interval Censored or Left Censored Data
    if (length(index)>0){
      gfw2$gamma[index]=dwg.fun(w2[index],q[index],y[index,2],x$gamma[index,],par$gamma,x$alpha[index,],par$alpha)
      gfw2$alpha[index]=dwa.fun(w2[index],q[index],y[index,2],x$gamma[index,],par$gamma,x$alpha[index,],par$alpha)
      if(!is.fix.par$q){gfw2$q[index]=dwq.fun(w2[index],q[index])}
      gSw2$gamma[index]=Swg.fun(w2[index],q[index],y[index,2],x$gamma[index,],par$gamma,x$alpha[index,],par$alpha)
      gSw2$alpha[index]=Swa.fun(w2[index],q[index],y[index,2],x$gamma[index,],par$gamma,x$alpha[index,],par$alpha)
      if(!is.fix.par$q){gSw2$q[index]=Swq.fun(w2[index],q[index])}
    }

    if(outtype=="HLLF"){
      ##### Second Derivative of Density and Survival Functions of w
      ### "hfw1","hSw1","hfw2","hSw2"
      hfw1=list(gamma=rep(NA,nobs),alpha=rep(NA,nobs),q=rep(NA,nobs),
                ga=rep(NA,nobs),gq=rep(NA,nobs),aq=rep(NA,nobs))
      hSw1=list(gamma=rep(NA,nobs),alpha=rep(NA,nobs),q=rep(NA,nobs),
                ga=rep(NA,nobs),gq=rep(NA,nobs),aq=rep(NA,nobs))
      hfw2=list(gamma=rep(NA,nobs),alpha=rep(NA,nobs),q=rep(NA,nobs),
                ga=rep(NA,nobs),gq=rep(NA,nobs),aq=rep(NA,nobs))
      hSw2=list(gamma=rep(NA,nobs),alpha=rep(NA,nobs),q=rep(NA,nobs),
                ga=rep(NA,nobs),gq=rep(NA,nobs),aq=rep(NA,nobs))

      index=which(y[3]==1|y[3]==3|y[3]==2) # Exact, Interval Censored or Right Censored Data
      if (length(index)>0){
        hfw1$gamma[index]=dwgg.fun(w1[index],q[index],y[index,1],x$gamma[index,],par$gamma,x$alpha[index,],par$alpha)
        hfw1$ga[index]=dwga.fun(w1[index],q[index],y[index,1],x$gamma[index,],par$gamma,x$alpha[index,],par$alpha)
        if(!is.fix.par$q) hfw1$gq[index]=dwgq.fun(w1[index],q[index],y[index,1],x$gamma[index,],par$gamma,x$alpha[index,],par$alpha)
        hfw1$alpha[index]=dwaa.fun(w1[index],q[index],y[index,1],x$gamma[index,],par$gamma,x$alpha[index,],par$alpha)
        if(!is.fix.par$q) hfw1$aq[index]=dwaq.fun(w1[index],q[index],y[index,1],x$gamma[index,],par$gamma,x$alpha[index,],par$alpha)
        if(!is.fix.par$q) hfw1$q[index]=dwqq.fun(w1[index],q[index])

        hSw1$gamma[index]=Swgg.fun(w1[index],q[index],y[index,1],x$gamma[index,],par$gamma,x$alpha[index,],par$alpha)
        hSw1$ga[index]=Swga.fun(w1[index],q[index],y[index,1],x$gamma[index,],par$gamma,x$alpha[index,],par$alpha)
        if(!is.fix.par$q) hSw1$gq[index]=Swgq.fun(w1[index],q[index],y[index,1],x$gamma[index,],par$gamma,x$alpha[index,],par$alpha)

        hSw1$alpha[index]=Swaa.fun(w1[index],q[index],y[index,1],x$gamma[index,],par$gamma,x$alpha[index,],par$alpha)
        if(!is.fix.par$q) hSw1$aq[index]=Swaq.fun(w1[index],q[index],y[index,1],x$gamma[index,],par$gamma,x$alpha[index,],par$alpha)
        if(!is.fix.par$q) hSw1$q[index]=Swqq.fun(w1[index],q[index])        
      }

      index=which(y[3]==3|y[3]==4) # Interval or Left Censored Data
      if (length(index)>0){
        hfw2$gamma[index]=dwgg.fun(w2[index],q[index],y[index,2],x$gamma[index,],par$gamma,x$alpha[index,],par$alpha)
        hfw2$ga[index]=dwga.fun(w2[index],q[index],y[index,2],x$gamma[index,],par$gamma,x$alpha[index,],par$alpha)
        if(!is.fix.par$q) hfw2$gq[index]=dwgq.fun(w2[index],q[index],y[index,2],x$gamma[index,],par$gamma,x$alpha[index,],par$alpha)
        hfw2$alpha[index]=dwaa.fun(w2[index],q[index],y[index,2],x$gamma[index,],par$gamma,x$alpha[index,],par$alpha)
        if(!is.fix.par$q) hfw2$aq[index]=dwaq.fun(w2[index],q[index],y[index,2],x$gamma[index,],par$gamma,x$alpha[index,],par$alpha)
        if(!is.fix.par$q) hfw2$q[index]=dwqq.fun(w2[index],q[index])
        hSw2$gamma[index]=Swgg.fun(w2[index],q[index],y[index,2],x$gamma[index,],par$gamma,x$alpha[index,],par$alpha)
        hSw2$ga[index]=Swga.fun(w2[index],q[index],y[index,2],x$gamma[index,],par$gamma,x$alpha[index,],par$alpha)
        if(!is.fix.par$q) hSw2$gq[index]=Swgq.fun(w2[index],q[index],y[index,2],x$gamma[index,],par$gamma,x$alpha[index,],par$alpha)
        hSw2$alpha[index]=Swaa.fun(w2[index],q[index],y[index,2],x$gamma[index,],par$gamma,x$alpha[index,],par$alpha)
        if(!is.fix.par$q) hSw2$aq[index]=Swaq.fun(w2[index],q[index],y[index,2],x$gamma[index,],par$gamma,x$alpha[index,],par$alpha)
        if(!is.fix.par$q) hSw2$q[index]=Swqq.fun(w2[index],q[index])
      }
    }
  }

  ########################################################################
  ##### Log-likelihood Contribution with Each Obervation with Each Subject
  ##### and It's First And Second Derivative
  ########################################################################
  if(outtype=="LLF"){
    ##### Log-likelihood Contribution with Each Obervation with Each Subject
    LLC=rep(0,nobs)

    ## Exact Data
    index=which(y[3]==1)
    if (length(index)>0) LLC[index]=log(p1[index])+log(fw1[index])-log(alpha[index])-y[index,1]
    ## Interval Censored Data
    index=which(y[3]==3)
    if (length(index)>0) LLC[index]=log(p1[index])+log(Sw1[index]-Sw2[index])
    ## Left Censored Data
    index=which(y[3]==4)
    if (length(index)>0) LLC[index]=log(p1[index]*(1-Sw2[index]))
    ## Right Censored Data
    index=which(y[3]==2)
    if (length(index)>0) LLC[index]=log(p1[index]*Sw1[index]+p2[index])

    ##### Weighting and Negative Taking of Log-Likelihood Contribution with Each Subject
    LLC=-LLC*weight

    ##### Summation of Log-Likelihood Contribution from All Subjects
    TLLC=sum(LLC)

    THGNlogLC=TLLC
  }else if(outtype=="GLLF"){
    ##### First Derivative Log-likelihood Contribution with Each Subject
    GLLC=list(beta=rep(0,nobs),gamma=rep(0,nobs),alpha=rep(0,nobs),q=rep(0,nobs))

    ## Exact Data
    index=which(y[3]==1)
    if (length(index)>0){
      base1=fw1[index]
      base1g=gfw1$gamma[index]
      base1a=gfw1$alpha[index]
      base1q=gfw1$q[index]
      if (is.reg$beta){ # beta
        GLLC$beta[index]=gp1$beta[index]/p1[index]
      }
      GLLC$gamma[index]=base1g/base1
      GLLC$alpha[index]=base1a/base1-Aa.fun(x$alpha,par$alpha)[index]/A.fun(x$alpha,par$alpha)[index]
      if(!is.fix.par$q) GLLC$q[index]=base1q/base1
    }
    ## Interval Censored Data
    index=which(y[3]==3)
    if (length(index)>0){
      base1=Sw1[index]-Sw2[index]
      base1g=gSw1$gamma[index]-gSw2$gamma[index]
      base1a=gSw1$alpha[index]-gSw2$alpha[index]      
      base1q=gSw1$q[index]-gSw2$q[index]
      if (is.reg$beta){ # beta
        GLLC$beta[index]=ip1[index]*gp1$beta[index]/p1[index]
      }
      GLLC$gamma[index]=base1g/base1
      GLLC$alpha[index]=base1a/base1
      if(!is.fix.par$q) GLLC$q[index]=base1q/base1
    }
    ## Left Censored Data
    index=which(y[3]==4)
    if (length(index)>0){
      base1=p1[index]*(1-Sw2[index])
      base1b=ip1[index]*gp1$beta[index]*(1-Sw2[index])
      base1g=p1[index]*(-gSw2$gamma[index])
      base1a=p1[index]*(-gSw2$alpha[index])
      base1q=p1[index]*(-gSw2$q[index])
      if (is.reg$beta){ # beta
        GLLC$beta[index]=base1b/base1
      }
      GLLC$gamma[index]=base1g/base1
      GLLC$alpha[index]=base1a/base1
      if(!is.fix.par$q) GLLC$q[index]=base1q/base1
    }
    ## Right Censored Data
    index=which(y[3]==2)
    if (length(index)>0){
      base1=p1[index]*Sw1[index]+p2[index]
      base1b=ip1[index]*gp1$beta[index]*Sw1[index]+ip2[index]*gp2$beta[index]
      base1g=p1[index]*gSw1$gamma[index]
      base1a=p1[index]*gSw1$alpha[index]
      base1q=p1[index]*gSw1$q[index]
      if (is.reg$beta){ # beta
        GLLC$beta[index]=base1b/base1
      }
      GLLC$gamma[index]=base1g/base1
      GLLC$alpha[index]=base1a/base1
      if(!is.fix.par$q) GLLC$q[index]=base1q/base1
    }

    ##### Weighting and Negative Taking of First Derivative Log-Likelihood Contribution with Each Subject
    if (is.reg$beta) GLLC$beta=-GLLC$beta*weight*x$beta
    GLLC$gamma=-GLLC$gamma*weight*x$gamma
    GLLC$alpha=-GLLC$alpha*weight*x$alpha
    if(!is.fix.par$q) GLLC$q=-GLLC$q*weight*x$q

    ##### Summation of First Derivative Log-Likelihood Contribution from All Subjects
    TGLLC=NULL
    if (is.reg$beta) TGLLC=c(TGLLC,as.vector(apply(GLLC$beta,2,sum))) # beta
    TGLLC=c(TGLLC,as.vector(apply(GLLC$gamma,2,sum))) # gamma
    TGLLC=c(TGLLC,as.vector(apply(GLLC$alpha,2,sum))) # alpha
    if(!is.fix.par$q) TGLLC=c(TGLLC,as.vector(apply(GLLC$q,2,sum))) # q

    THGNlogLC=TGLLC
  }else if(outtype=="HLLF"){
    ##### Second Derivative Log-likelihood Contribution with Each Subject
    HLLC=list(beta=rep(0,nobs),gamma=rep(0,nobs),alpha=rep(0,nobs),q=rep(0,nobs),
              bg=rep(0,nobs),ba=rep(0,nobs),bq=rep(0,nobs),
              ga=rep(0,nobs),gq=rep(0,nobs),aq=rep(0,nobs))

    ## Exact Data
    index=which(y[3]==1)
    if (length(index)>0){
      base1=fw1[index]
      base1g=gfw1$gamma[index]
      base1a=gfw1$alpha[index]
      base1q=gfw1$q[index]                  
      base1gg=hfw1$gamma[index]
      base1ga=hfw1$ga[index]
      base1gq=hfw1$gq[index]
      base1aa=hfw1$alpha[index]
      base1aq=hfw1$aq[index]
      base1qq=hfw1$q[index]
      if (is.reg$beta) HLLC$beta[index]=(hp1$beta[index]*p1[index]-gp1$beta[index]^2)/p1[index]^2
      HLLC$gamma[index]=(base1gg*base1-base1g^2)/base1^2
      HLLC$ga[index]=(base1ga*base1-base1g*base1a)/base1^2
      if(!is.fix.par$q) HLLC$gq[index]=(base1gq*base1-base1g*base1q)/base1^2
      HLLC$alpha[index]=(base1aa*base1-base1a^2)/base1^2-
                        (Aaa.fun(x$alpha,par$alpha)[index]*A.fun(x$alpha,par$alpha)[index]-Aa.fun(x$alpha,par$alpha)[index]^2)/A.fun(x$alpha,par$alpha)[index]^2
      if(!is.fix.par$q) HLLC$aq[index]=(base1aq*base1-base1a*base1q)/base1^2
      if(!is.fix.par$q) HLLC$q[index]=(base1qq*base1-base1q^2)/base1^2
    }
    ## Right Censored Data
    index=which(y[3]==2)
    if (length(index)>0){
      base1=p1[index]*Sw1[index]+p2[index]
      base1b=ip1[index]*gp1$beta[index]*Sw1[index]+ip2[index]*gp2$beta[index]
      base1g=p1[index]*gSw1$gamma[index]
      base1a=p1[index]*gSw1$alpha[index]
      base1q=p1[index]*gSw1$q[index]
      base1bb=ip1[index]*hp1$beta[index]*Sw1[index]+ip2[index]*hp2$beta[index]
      base1bg=ip1[index]*gp1$beta[index]*gSw1$gamma[index]
      base1ba=ip1[index]*gp1$beta[index]*gSw1$alpha[index]
      base1bq=ip1[index]*gp1$beta[index]*gSw1$q[index]
      base1gg=p1[index]*hSw1$gamma[index]
      base1ga=p1[index]*hSw1$ga[index]
      base1gq=p1[index]*hSw1$gq[index]
      base1aa=p1[index]*hSw1$alpha[index]
      base1aq=p1[index]*hSw1$aq[index]
      base1qq=p1[index]*hSw1$q[index]
      if (is.reg$beta){
        HLLC$beta[index]=(base1bb*base1-base1b^2)/base1^2
        HLLC$bg[index]=(base1bg*base1-base1b*base1g)/base1^2
        HLLC$ba[index]=(base1ba*base1-base1b*base1a)/base1^2
        if(!is.fix.par$q) HLLC$bq[index]=(base1bq*base1-base1b*base1q)/base1^2
      }
      HLLC$gamma[index]=(base1gg*base1-base1g^2)/base1^2
      HLLC$ga[index]=(base1ga*base1-base1g*base1a)/base1^2
      if(!is.fix.par$q) HLLC$gq[index]=(base1gq*base1-base1g*base1q)/base1^2
      HLLC$alpha[index]=(base1aa*base1-base1a^2)/base1^2
      if(!is.fix.par$q) HLLC$aq[index]=(base1aq*base1-base1a*base1q)/base1^2
      if(!is.fix.par$q) HLLC$q[index]=(base1qq*base1-base1q^2)/base1^2
    }
    ## Interval Censored Data
    index=which(y[3]==3)
    if (length(index)>0){
      base1=Sw1[index]-Sw2[index]
      base1g=gSw1$gamma[index]-gSw2$gamma[index]
      base1a=gSw1$alpha[index]-gSw2$alpha[index]      
      base1q=gSw1$q[index]-gSw2$q[index]
      base1gg=hSw1$gamma[index]-hSw2$gamma[index]
      base1ga=hSw1$ga[index]-hSw2$ga[index]
      base1gq=hSw1$gq[index]-hSw2$gq[index]
      base1aa=hSw1$alpha[index]-hSw2$alpha[index]
      base1aq=hSw1$aq[index]-hSw2$aq[index]
      base1qq=hSw1$q[index]-hSw2$q[index]
      if (is.reg$beta){
        HLLC$beta[index]=(hp1$beta[index]*p1[index]-gp1$beta[index]^2)/p1[index]^2
      }
      HLLC$gamma[index]=(base1gg*base1-base1g^2)/base1^2
      HLLC$ga[index]=(base1ga*base1-base1g*base1a)/base1^2
      if(!is.fix.par$q) HLLC$gq[index]=(base1gq*base1-base1g*base1q)/base1^2
      HLLC$alpha[index]=(base1aa*base1-base1a^2)/base1^2
      if(!is.fix.par$q) HLLC$aq[index]=(base1aq*base1-base1a*base1q)/base1^2
      if(!is.fix.par$q) HLLC$q[index]=(base1qq*base1-base1q^2)/base1^2
    }
    ## Left Censored Data
    index=which(y[3]==4)
    if (length(index)>0){
      base1=p1[index]*(1-Sw2[index])
      base1b=ip1[index]*gp1$beta[index]*(1-Sw2[index])
      base1g=p1[index]*(-gSw2$gamma[index])
      base1a=p1[index]*(-gSw2$alpha[index])
      base1q=p1[index]*(-gSw2$q[index])
      base1bb=ip1[index]*hp1$beta[index]*(1-Sw2[index])
      base1bg=ip1[index]*gp1$beta[index]*(-gSw2$gamma[index])
      base1ba=ip1[index]*gp1$beta[index]*(-gSw2$alpha[index])
      base1bq=ip1[index]*gp1$beta[index]*(-gSw2$q[index])
      base1gg=p1[index]*(-hSw2$gamma[index])
      base1ga=p1[index]*(-hSw2$ga[index])
      base1gq=p1[index]*(-hSw2$gq[index])
      base1aa=p1[index]*(-hSw2$alpha[index])
      base1aq=p1[index]*(-hSw2$aq[index])
      base1qq=p1[index]*(-hSw2$q[index])
      if (is.reg$beta){
        HLLC$beta[index]=(base1bb*base1-base1b^2)/base1^2
        HLLC$bg[index]=(base1bg*base1-base1b*base1g)/base1^2
        HLLC$ba[index]=(base1ba*base1-base1b*base1a)/base1^2
        if(!is.fix.par$q) HLLC$bq[index]=(base1bq*base1-base1b*base1q)/base1^2
      }
      HLLC$gamma[index]=(base1gg*base1-base1g^2)/base1^2
      HLLC$ga[index]=(base1ga*base1-base1g*base1a)/base1^2
      if(!is.fix.par$q) HLLC$gq[index]=(base1gq*base1-base1g*base1q)/base1^2
      HLLC$alpha[index]=(base1aa*base1-base1a^2)/base1^2
      if(!is.fix.par$q) HLLC$aq[index]=(base1aq*base1-base1a*base1q)/base1^2
      if(!is.fix.par$q) HLLC$q[index]=(base1qq*base1-base1q^2)/base1^2
    }

    ##### Weighting and Negative Taking of Second Derivative Log-Likelihood Contribution with Each Subject
    if (is.reg$beta){
      HLLC$beta=-HLLC$beta*weight*CMultiply(x$beta)
      HLLC$bg=-HLLC$bg*weight*CMultiply(x$beta,x$gamma)
      HLLC$ba=-HLLC$ba*weight*CMultiply(x$beta,x$alpha)
      if(!is.fix.par$q) HLLC$bq=-HLLC$bq*weight*CMultiply(x$beta,x$q)
    }
    HLLC$gamma=-HLLC$gamma*weight*CMultiply(x$gamma)
    HLLC$ga=-HLLC$ga*weight*CMultiply(x$gamma,x$alpha)
    if(!is.fix.par$q) HLLC$gq=-HLLC$gq*weight*CMultiply(x$gamma,x$q)
    HLLC$alpha=-HLLC$alpha*weight*CMultiply(x$alpha)
    if(!is.fix.par$q) HLLC$aq=-HLLC$aq*weight*CMultiply(x$alpha,x$q)
    if(!is.fix.par$q) HLLC$q=-HLLC$q*weight*CMultiply(x$q)

    ##### Summation of Second Derivative Log-Likelihood Contribution from All Subjects
    loc2=c(length(par$beta),length(par$beta)+length(par$gamma),length(par$beta)+length(par$gamma)+length(par$alpha),
           length(par$beta)+length(par$gamma)+length(par$alpha)+length(par$q))
    loc1=loc2-c(length(par$beta),length(par$gamma),length(par$alpha),length(par$q))+1
    dim=0
    if (is.reg$beta) dim=dim+length(par$beta)
    if (is.reg$gamma) dim=dim+length(par$gamma)
    if (is.reg$alpha) dim=dim+length(par$alpha)
    if (is.reg$q){
      if(!is.fix.par$q) dim=dim+length(par$q)
    }
    THLLC=matrix(0,dim,dim)
    ### beta
    if (is.reg$beta){
      # beta
      buffer=as.vector(apply(HLLC$beta,2,sum));k=1
      for (i in loc1[1]:loc2[1]){
        for(j in i:loc2[1]){
          THLLC[i,j]=buffer[k];k=k+1
        }
      }
      # bg
      buffer=as.vector(apply(HLLC$bg,2,sum));k=1
      for (i in loc1[1]:loc2[1]){
        for(j in loc1[2]:loc2[2]){
          THLLC[i,j]=buffer[k];k=k+1
        }
      }
      # ba
      buffer=as.vector(apply(HLLC$ba,2,sum));k=1
      for (i in loc1[1]:loc2[1]){
        for(j in loc1[3]:loc2[3]){
          THLLC[i,j]=buffer[k];k=k+1
        }
      }
      # bq
      if(!is.fix.par$q){ 
        buffer=as.vector(apply(HLLC$bq,2,sum));k=1
        for (i in loc1[1]:loc2[1]){
          for(j in loc1[4]:loc2[4]){
            THLLC[i,j]=buffer[k];k=k+1
          }
        }
      }
    }

    ### gamma
    buffer=as.vector(apply(HLLC$gamma,2,sum));k=1
    for (i in loc1[2]:loc2[2]){ # gamma
      for(j in i:loc2[2]){
        THLLC[i,j]=buffer[k];k=k+1
      }
    }
    # ga
    buffer=as.vector(apply(HLLC$ga,2,sum));k=1
    for (i in loc1[2]:loc2[2]){ 
      for(j in loc1[3]:loc2[3]){
        THLLC[i,j]=buffer[k];k=k+1
      }
    }
    # gq
    if(!is.fix.par$q){ 
      buffer=as.vector(apply(HLLC$gq,2,sum));k=1
      for (i in loc1[2]:loc2[2]){ 
        for(j in loc1[4]:loc2[4]){
          THLLC[i,j]=buffer[k];k=k+1
        }
      }
    }

    ### alpha
    buffer=as.vector(apply(HLLC$alpha,2,sum));k=1
    for (i in loc1[3]:loc2[3]){ # alpha
      for(j in i:loc2[3]){
        THLLC[i,j]=buffer[k];k=k+1
      }
    }
    if(!is.fix.par$q){ # aq
      buffer=as.vector(apply(HLLC$aq,2,sum));k=1
      for (i in loc1[3]:loc2[3]){ 
        for(j in loc1[4]:loc2[4]){
          THLLC[i,j]=buffer[k];k=k+1
        }
      }
    }

    ### q
    if(!is.fix.par$q){
      buffer=as.vector(apply(HLLC$q,2,sum));k=1
      for (i in loc1[4]:loc2[4]){ # q
        for(j in i:loc2[4]){THLLC[i,j]=buffer[k];k=k+1}}
    }
    THGNlogLC=t(THLLC)+THLLC;diag(THGNlogLC)=diag(THGNlogLC)/2
  }
  return(THGNlogLC)
}

#########################
### function: CMultiply()
#########################
CMultiply=function(A,B=NULL){
  A=as.matrix(A)
  if (is.null(B)){B=A;is.same=T
  }else{B=as.matrix(B);is.same=F}
  dimA=dim(A);dimB=dim(B)
  result=NULL
  if (is.same){
    for (i in 1:dimA[2]) result=cbind(result,A[,i]*B[,i:dimB[2]])
  }else{
    if (dimA[1]==dimB[1]){
      for (i in 1:dimA[2]) result=cbind(result,A[,i]*B)
    } 
  }
  return(result)
}
# CMultiply(x$alpha)
# CMultiply(x$alpha,x$beta)
