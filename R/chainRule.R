####################################################################################################
##### The functions used related to the log-likelihood function
####################################################################################################
############################
### Probability of binormial
############################
### "p1.fun"
p1.fun=function(x,par,is.reg){
  if(!is.reg$beta) return(0)
  reg=as.matrix(x$beta)%*%par$beta
  result=as.vector(1/(1+exp(-reg)))
  return(result)
}

##################
### Scale function
##################
### "A.fun" (Scale Component)
A.fun=function(x.alpha,alpha){
  result=as.vector(exp(as.matrix(x.alpha)%*%as.vector(alpha)))
  return(result)
}

##################
### Shape function
##################
### "Q.fun" (Shape Component)
Q.fun=function(x.q,q){
  if(q[1]=="logistic") return(rep("logistic",dim(x.q)[1]))
  result=as.vector(as.matrix(x.q)%*%as.vector(q))
  return(result)
}

##############
### w function
##############
### "w.fun"
w.fun=function(y,x.gamma,gamma,x.alpha,alpha){
  xgamma=as.matrix(x.gamma)%*%gamma;xalpha=as.matrix(x.alpha)%*%alpha
  result=as.vector((y-xgamma)/exp(xalpha))
  if(length(y)==1) y=rep(y,length(result))
  index=which(abs(y)==Inf);if(length(index)>0) result[index]=y[index]
  return(result)
}

####################################################################
### Density, survival and quantile function of logistic distribution
####################################################################
### "dlog","slog","qlog"
dlog=function(w){
  result=rep(NA,length(w))
  index=which(w<=0);if(length(index)>0) result[index]=exp(w[index])/(1+exp(w[index]))^2
  index=which(w>0);if(length(index)>0) result[index]=1/(exp(-w[index])+2+exp(w[index]))
  index=is.na(result);if(length(index)>0) result[index]=0
  return(result)
}
slog=function(w){
  result=1/(1+exp(w))
  return(result)
}
qlog=function(prob){
  result=rep(NA,length(prob))
  index=which(prob>=0 & prob<=1);if(length(index)>0) result[index]=log(prob[index]/(1-prob[index]))
  return(result)
}

##################################################################
### Density, survival and quantile function of normal distribution
##################################################################
### "dnorm","snorm","qnorm"
#dnorm=>default function
snorm=function(w){
  result=1-pnorm(w)
  return(result)
}
#qnorm=>default function

#################################################################################
### Density, survival and quantile function of generalized log-gamma distribution
#################################################################################
### "dglgd","sglgd","qglgd"
dglgd=function(w,q){
  if(length(w)!=length(q)){
    lenw=length(w);lenq=length(q);len=max(lenw,lenq)
    w=rep(w,ceiling(len/lenw))[1:len];q=rep(q,ceiling(len/lenq))[1:len]
  }
  q2=q^(-2)
  result=exp(-lgamma(q2)+log(abs(q))+q2*log(q2)+q2*(q*w-exp(q*w)))
  index=which(q==0);if(length(index)>0) result[index]=dnorm(w[index])
  index=is.na(result);if(length(index)>0) result[index]=0
  return(result)
}
if(T){
  sglgd=function(w,q){
    if(length(w)!=length(q)){
      lenw=length(w);lenq=length(q);len=max(lenw,lenq)
      w=rep(w,ceiling(len/lenw))[1:len];q=rep(q,ceiling(len/lenq))[1:len]
    }
    ev=exp(q*w-log(q^2));q2=q^(-2)
    result=IGamma(x=ev,p=q2)[,6]
    index=which(q>0);if(length(index)>0) result[index]=1-result[index]
    index=which(q==0);if(length(index)>0) result[index]=1-pnorm(w[index])
    index=which(w>0 & result==1);if(length(index)>0) result[index]=0 # Adjust for q>0
    index=which(w<0 & result==0);if(length(index)>0) result[index]=1 # Adjust for q<0
    return(result)
  }
}else{
  sglgd=function(w,q){
    if(length(w)!=length(q)){
      lenw=length(w);lenq=length(q);len=max(lenw,lenq)
      w=rep(w,ceiling(len/lenw))[1:len];q=rep(q,ceiling(len/lenq))[1:len]
    }
    ev=exp(q*w-log(q^2));q2=q^(-2)
    result=rep(NA,length(w))
    index=which(q!=0);if(length(index)>0) result[index]=pgamma(ev[index],shape=q2[index],scale=1)
    index=which(q>0);if(length(index)>0) result[index]=1-result[index]
    index=which(q==0);if(length(index)>0) result[index]=1-pnorm(w[index])
    return(result)
  }
}
qglgd=function(prob,q){
  if(length(prob)!=length(q)){
    lenprob=length(prob);lenq=length(q);len=max(lenprob,lenq)
    prob=rep(prob,ceiling(len/lenprob))[1:len];q=rep(q,ceiling(len/lenq))[1:len]
  }
  q2=q^(-2)
  index=which(q<0);if(length(index)>0) prob[index]=1-prob[index]
  result=rep(NA,length(prob))
  index=which(q==0);if(length(index)>0) result[index]=qnorm(p=prob[index])
  index=which(q!=0);if(length(index)>0) result[index]=(log(qgamma(p=prob[index],shape=q2[index],scale=1))+log(q[index]^2))/q[index]
  return(result)
}

#####################################################################
### Density, survival and quantile function of w (error distribution)
#####################################################################
### "dw.fun", "Sw.fun", "qw.fun"
dw.fun=function(w,q){
  if(length(w)!=length(q)){
    lenw=length(w);lenq=length(q);len=max(lenw,lenq)
    w=rep(w,ceiling(len/lenw))[1:len];q=rep(q,ceiling(len/lenq))[1:len]
  }
  result=rep(NA,length(w))
  # logistic distribution
  index=which(q=="logistic");if(length(index)>0) result[index]=dlog(w[index])
  # normal distribution
  if(is.numeric(q)){
    index=which(q==0);if(length(index)>0) result[index]=dnorm(w[index])
    index=which(abs(q)<1e-05);if(length(index)>0) result[index]=dnorm(w[index])
  }
  # generalized log-gamma distribution
  if(is.numeric(q)){
    index=which(q!=0 & abs(q)>=1e-05);if(length(index)>0) result[index]=dglgd(w[index],q[index])
  }
  return(result)
}
Sw.fun=function(w,q){
  if(length(w)!=length(q)){
    lenw=length(w);lenq=length(q);len=max(lenw,lenq)
    w=rep(w,ceiling(len/lenw))[1:len];q=rep(q,ceiling(len/lenq))[1:len]
  }
  result=rep(NA,length(q))
  # logistic distribution
  index=which(q=="logistic");if(length(index)>0) result[index]=slog(w[index])
  # normal distribution
  if(is.numeric(q)){
    index=which(q==0);if(length(index)>0) result[index]=1-pnorm(w[index])
    index=which(abs(q)<1e-5);if(length(index)>0) result[index]=1-pnorm(w[index])
  }
  # generalized log-gamma distribution
  if(is.numeric(q)){
    index=which(q!=0 & abs(q)>=1e-5);if(length(index)>0) result[index]=sglgd(w[index],q[index])
  }
  return(result)
}
qw.fun=function(prob,q){
  if(length(prob)!=length(q)){
    lenprob=length(prob);lenq=length(q);len=max(lenprob,lenq)
    prob=rep(prob,ceiling(len/lenprob))[1:len];q=rep(q,ceiling(len/lenq))[1:len]
  }
  result=rep(NA,length(prob))
  # logistic distribution
  index=which(q=="logistic");if(length(index)>0) result[index]=qlog(prob[index])
  # normal distribution
  if(is.numeric(q)){
    index=which(q==0);if(length(index)>0) result[index]=qnorm(prob[index])
  }
  # generalized log-gamma distribution
  if(is.numeric(q)){
    index=which(q!="logistic" & q!=0);if(length(index)>0) result[index]=qglgd(prob[index],as.numeric(q[index]))
  }
  return(result)
}

###########################################################
### Density and Survival Function for T (time distribution)
###########################################################
### "ft", "St"
ft=function(time,xgamma,xalpha,xq){
  if(length(time)!=length(xq)) xq=rep(xq[1],length(time))
  w=(log(time)-xgamma)/exp(xalpha)
  result=dw.fun(w,xq)/(exp(xalpha)*time)
  index=which(time==0)
  if(length(index)>0) result[index]=0
  return(result)
}
St=function(time,xgamma,xalpha,xq){
  if(length(time)!=length(xq))xq=rep(xq[1],length(time))
  w=(log(time)-xgamma)/exp(xalpha)
  result=Sw.fun(w,xq)
  return(result)
}
############
### Example:
############
if(F){
 xgamma=3.8;xalpha=-1.049822;xq=1
 ft(0:40,xgamma=xgamma,xalpha=xalpha,xq=xq)
 St(0:40,xgamma=xgamma,xalpha=xalpha,xq=xq)
 integrate(ft,0,40,xgamma=xgamma,xalpha=xalpha,xq=xq)
 1-St(40,xgamma=xgamma,xalpha=xalpha,xq=xq)
}

############################################################################
# Incomplete gamma integral and it's derivative through normal approximation
# IGammaNormal(x,p)[1] : 1st derivative of IGammaNormal(x,p) w.r.t. x
# IGammaNormal(x,p)[2] : 2nd derivative of IGammaNormal(x,p) w.r.t. x^2
# IGammaNormal(x,p)[3] : 1st derivative of IGammaNormal(x,p) w.r.t. p
# IGammaNormal(x,p)[4] : 2nd derivative of IGammaNormal(x,p) w.r.t. p^2
# IGammaNormal(x,p)[5] : 2nd derivative of IGammaNormal(x,p) w.r.t. x and p
# IGammaNormal(x,p)[6] : value of IGammaNormal(x,p)
############################################################################
IGammaNormal=function(x,p){
  num=length(x)
  result=matrix(NA,num,6)
  y=3*p^(1/6)*x^(1/3)+(1/3)*p^(-1/2)-3*p^(1/2)
  y=3*sqrt(p)*((x/p)^(1/3)+1/(9*p)-1)
  dyx=p^(1/6)*x^(-2/3)
  dyxx=(-2/3)*p^(1/6)*x^(-5/3)
  dyp=(1/2)*p^(-5/6)*x^(1/3)-(1/6)*p^(-3/2)-(3/2)*p^(-1/2)
  dypp=(-5/12)*p^(-11/6)*x^(1/3)+(1/4)*p^(-5/2)+(3/4)*p^(-3/2)
  dyxp=(1/6)*p^(-5/6)*x^(-2/3)
  result[,1]=dnorm(y)*dyx
  result[,2]=(-y*dnorm(y))*dyx^2+dnorm(y)*dyxx
  result[,3]=dnorm(y)*dyp
  result[,4]=(-y*dnorm(y))*dyp^2+dnorm(y)*dypp
  result[,5]=(-y*dnorm(y))*dyx*dyp+dnorm(y)*dyxp
  result[,6]=pnorm(y)
  return(result)
}
### Example:
#IGammaNormal(x=c(1,2),p=c(15,15))
#x=c(1,2);p=c(15,15);IGammaNormal(x=c(30),p=c(50))
###############################################################
# Incomplete gamma integral and it's derivative
# IGamma(x,p)[1] : 1st derivative of IGamma(x,p) w.r.t. x
# IGamma(x,p)[2] : 2nd derivative of IGamma(x,p) w.r.t. x^2
# IGamma(x,p)[3] : 1st derivative of IGamma(x,p) w.r.t. p
# IGamma(x,p)[4] : 2nd derivative of IGamma(x,p) w.r.t. p^2
# IGamma(x,p)[5] : 2nd derivative of IGamma(x,p) w.r.t. x and p
# IGamma(x,p)[6] : value of IGamma(x,p)
###############################################################
IGamma=function(x,p,plimit=10000000000){
  index=which(x>1e+99);if(length(index)>0)x[index]=1e+99
  index=which(p>1e+99);if(length(index)>0)p[index]=1e+99
  num=length(x);d=rep(0,num*6);plimit=plimit;ifault=0
  result=.Fortran("digamiv",as.integer(num),as.double(d),as.double(x),as.double(p),as.double(plimit),as.integer(ifault))
  result=matrix(result[[2]],length(x),6,byrow=T)
  index=which(p>plimit);if(length(index)>0) result[index,]=IGammaNormal(x[index],p[index])
  return(result)
}
### Example:
#IGamma(x=c(1,2,1,10000,1000),p=c(1,2,10000,15,1000))

####################################################################################################
##### The functions used related to the 1st derivative of the log-likelihood function
####################################################################################################
#########################################################
### 1st derivative of the probability of nested-binormial
#########################################################
### "p1b.fun"
p1b.fun=function(x,par,is.reg){
  if(!is.reg$beta) return(0)
  p=p1.fun(x,par,is.reg)
  result=p*(1-p)
  return(result)
}
########################################
### 1st derivative of the Scale function
########################################
### "Aa.fun"
Aa.fun=function(x.alpha,alpha){
  result=A.fun(x.alpha,alpha)
  return(result)
}
################################
### 1st derivative of w function
################################
### "wg.fun","wa.fun"
wg.fun=function(y,x.gamma,gamma,x.alpha,alpha){
  wfun=w.fun(y,x.gamma,gamma,x.alpha,alpha)
  xalpha=as.matrix(x.alpha)%*%alpha
  result=as.vector(-1/exp(xalpha))
  index=which(abs(wfun)==Inf);if(length(index)>0) result[index]=0
  return(result)
}
wa.fun=function(y,x.gamma,gamma,x.alpha,alpha){
  wfun=w.fun(y,x.gamma,gamma,x.alpha,alpha)
  index=which(abs(wfun)==Inf);if(length(index)>0) wfun[index]=0
  result=as.vector(-wfun)
  return(result)
}

#############################################################################
### 1st derivative of density and survival function for logistic distribution
#############################################################################
### "dlogw","slogw"
dlogw=function(w){
  result=rep(NA,length(w))
  index=which(w<=0);if (length(index)>0) result[index]=((1-exp(w[index]))/(1+exp(w[index])))*dlog(w[index])
  index=which(w>0);if (length(index)>0) result[index]=((exp(-w[index])-1)/(exp(-w[index])+1))*dlog(w[index])
  return(result)
}
slogw=function(w){
  result=-dlog(w)
  return(result)
}
###########################################################################
### 1st derivative of density and survival function for normal distribution
###########################################################################
### "dnormw","snormw"
dnormw=function(w){
  result=-w*dnorm(w)
  index=is.na(result);if(length(index)>0) result[index]=0
  return(result)
}
snormw=function(w){
  result=-dnorm(w)
  return(result)
}
##########################################################################################
### 1st derivative of density and survival function for generalized log-gamma distribution
##########################################################################################
### "dglgdw","dglgdq","sglgdw","sglgdq"
dglgdw=function(w,q){
  result=q^(-1)*(1-exp(q*w))*dglgd(w,q)
  index=is.na(result);if(length(index)>0) result[index]=0
  return(result)
}
dglgdq=function(w,q){
  dig=digamma(1/q^2);dglgd=dglgd(w,q)
  result=as.vector((q^(-3))*(-2+2*exp(q*w)+q^2-q*w-q*w*exp(q*w)+2*log(q^2)+2*dig)*dglgd)
  index=is.na(result);if(length(index)>0) result[index]=0
  return(result)
}
sglgdw=function(w,q){
  result=-dglgd(w,q)
  return(result)
}
sglgdq=function(w,q){
  sign.q=q
  index=which(q>0);if(length(index)>0) sign.q[index]=1
  index=which(q<0);if(length(index)>0) sign.q[index]=-1
  ev=exp(q*w-log(q^2));q2=q^(-2)
  igamma=IGamma(x=ev,p=q2)
  result=as.vector(sign.q*((2/q^3)*igamma[,3]-igamma[,1]*ev*(w-2/q)))
  index=is.na(result);if(length(index)>0) result[index]=0
  return(result)
}

#####################################################################
### 1st derivative of density and survival for w (error distribution)
#####################################################################
### "dww.fun","dwg.fun","dwa.fun","dwq.fun"
### "Sww.fun","Swg.fun","Swa.fun","Swq.fun"
dww.fun=function(w,q){
  if(length(w)!=length(q)){
    lenw=length(w);lenq=length(q);len=max(lenw,lenq)
    w=rep(w,ceiling(len/lenw))[1:len];q=rep(q,ceiling(len/lenq))[1:len]
  }
  result=rep(NA,length(w))
  # logistic distribution
  index=which(q=="logistic");if(length(index)>0) result[index]=dlogw(w[index])
  # normal distribution
  if(is.numeric(q)){
    index=which(q==0);if(length(index)>0) result[index]=dnormw(w[index])
    index=which(abs(q)<1e-05);if(length(index)>0) result[index]=dnormw(w[index])
  }
  # generalized log-gamma distribution
  if(is.numeric(q)){
    index=which(q!=0 & abs(q)>=1e-05);if(length(index)>0) result[index]=dglgdw(w[index],q[index])
  }
  return(result)
}
dwg.fun=function(w,q,y,x.gamma,gamma,x.alpha,alpha){
  result=dww.fun(w,q)*wg.fun(y,x.gamma,gamma,x.alpha,alpha)
  return(result)
}
dwa.fun=function(w,q,y,x.gamma,gamma,x.alpha,alpha){
  result=dww.fun(w,q)*wa.fun(y,x.gamma,gamma,x.alpha,alpha)
  return(result)
}
dwq.fun=function(w,q){
  result=dglgdq(w,q)
  return(result)
}
###
Sww.fun=function(w,q){
  result=-dw.fun(w,q)
  return(result)
}
Swg.fun=function(w,q,y,x.gamma,gamma,x.alpha,alpha){
  result=Sww.fun(w,q)*wg.fun(y,x.gamma,gamma,x.alpha,alpha)
  return(result)
}
Swa.fun=function(w,q,y,x.gamma,gamma,x.alpha,alpha){
  result=Sww.fun(w,q)*wa.fun(y,x.gamma,gamma,x.alpha,alpha)
  return(result)
}
Swq.fun=function(w,q){
  result=sglgdq(w,q)
  return(result)
}

####################################################################################################
##### The functions used related to the 2st derivative of the log-likelihood function
####################################################################################################
#########################################################
### 2nd derivative of the probability of nested-binormial
#########################################################
### "p1bb.fun"
p1bb.fun=function(x,par,is.reg){
  if(!is.reg$beta) return(0)
  p=p1.fun(x,par,is.reg)
  result=p*(1-p)*(1-2*p)
  return(result)
}
########################################
### 2nd derivative of the Scale function
########################################
### "Aaa.fun"
Aaa.fun=function(x.alpha,alpha){
  result=A.fun(x.alpha,alpha)
  return(result)
}
################################
### 2nd derivative of w function
################################
### "wgg.fun","wga.fun","waa.fun"
wgg.fun=function(y,x.gamma,gamma,x.alpha,alpha){
  xalpha=as.vector(as.matrix(x.alpha)%*%alpha);xalpha[]=0
  result=xalpha
  return(result)
}
wga.fun=function(y,x.gamma,gamma,x.alpha,alpha){
  xalpha=as.matrix(x.alpha)%*%alpha
  result=as.vector(1/exp(xalpha))
  return(result)
}
waa.fun=function(y,x.gamma,gamma,x.alpha,alpha){
  wfun=w.fun(y,x.gamma,gamma,x.alpha,alpha)
  index=which(abs(wfun)==Inf);if(length(index)>0) wfun[index]=0
  result=as.vector(wfun)
  return(result)
}

#############################################################################
### 2nd derivative of density and survival function for logistic distribution
#############################################################################
### "dlogww","slogww"
dlogww=function(w){
  result=rep(NA,length(w))
  index=which(w<=0);if (length(index)>0) result[index]=(exp(-w[index])+exp(w[index])-4)*dlog(w[index])*dlog(w[index])
  index=which(w>0);if (length(index)>0) result[index]=(exp(-w[index])+exp(w[index])-4)*dlog(w[index])*dlog(w[index])
  return(result)
}
slogww=function(w){
  result=-dlogw(w)
  return(result)
}
###########################################################################
### 2nd derivative of density and survival function for normal distribution
###########################################################################
### "dnormw","snormww"
dnormww=function(w){
  result=(w*w-1)*dnorm(w)
  index=is.na(result);if(length(index)>0) result[index]=0
  return(result)
}
snormww=function(w){
  result=-dnormw(w)
  return(result)
}
##########################################################################################
### 2nd derivative of density and survival function for generalized log-gamma distribution
##########################################################################################
### "dglgdww","dglgdqq","dglgdwq"
### "sglgdww","sglgdqq","sglgdwq"
dglgdww=function(w,q){
  result=(q^(-2)*(1-exp(q*w))^2-exp(q*w))*dglgd(w,q)
  index=is.na(result);if(length(index)>0) result[index]=0
  return(result)
}
dglgdqq=function(w,q){
  dig=digamma(1/q^2);tri=trigamma(1/q^2)
  dglgd=dglgd(w,q);dglgdq=dglgdq(w,q)
  result=as.vector(dglgd*((dglgdq/dglgd)^2+
         (-q^(-2)-exp(q*w)*w^2*q^(-2)+2*q^(-3)*(w+2*w*exp(q*w)+2*q^(-1)-2*tri*q^(-3))+
         6*q^(-4)*(1-exp(q*w)-log(q^2)-dig))))
  index=is.na(result);if(length(index)>0) result[index]=0
  return(result)
}
dglgdwq=function(w,q){
  dglgd=dglgd(w,q);dglgdq=dglgdq(w,q)
  result=as.vector((1-exp(q*w))*q^(-1)*dglgdq+(exp(q*w)-1-q*w*exp(q*w))*q^(-2)*dglgd)
  index=is.na(result);if(length(index)>0) result[index]=0
  return(result)
}
###
sglgdww=function(w,q){
  result=-dglgdw(w,q)
  return(result)
}
sglgdqq=function(w,q){
  sign.q=q
  index=which(q>0);if(length(index)>0) sign.q[index]=1
  index=which(q<0);if(length(index)>0) sign.q[index]=-1
  dig=digamma(1/q^2);v=q*w-log(q^2)
  igamma=IGamma(x=exp(v),p=q^(-2))
  result=as.vector(sign.q*(-6*q^(-4)*igamma[,3]-(2*q^(-3))^2*igamma[,4]-igamma[,1]*exp(v)*q^(-2)*(
         (q^(-2)-exp(v))*(q*w-2)^2-4*q^(-2)*(v-dig)*(q*w-2)+2)))
  index=is.na(result);if(length(index)>0) result[index]=0
  return(result)
}
sglgdwq=function(w,q){
  result=-dglgdq(w,q)
  return(result)
}
#####################################################################
### 2nd derivative of density and survival for w (error distribution)
#####################################################################
### "dwww.fun","dwwq.fun","dwgg.fun","dwga.fun","dwgq.fun","dwaa.fun","dwaq.fun","dwqq.fun"
### "Swww.fun","Swwq.fun","Swgg.fun","Swga.fun","Swgq.fun","Swaa.fun","Swaq.fun","Swqq.fun"
dwww.fun=function(w,q){
  result=rep(NA,length(w))
  # log-logistic distribution
  index=which(q=="logistic");if(length(index)>0) result[index]=dlogww(w[index])
  # normal distribution
  if(is.numeric(q)){
    index=which(q==0);if(length(index)>0) result[index]=dnormww(w[index])
    index=which(abs(q)<1e-3);if(length(index)>0) result[index]=dnormww(w[index])
  }
  # generalized log-gamma distribution
  if(is.numeric(q)){
    index=which(q!=0 & abs(q)>1e-3);if(length(index)>0) result[index]=dglgdww(w[index],q[index])
  }
  return(result)
}
dwwq.fun=function(w,q){
  result=dglgdwq(w,q)
  return(result)
}
dwgg.fun=function(w,q,y,x.gamma,gamma,x.alpha,alpha){
  result=dwww.fun(w,q)*wg.fun(y,x.gamma,gamma,x.alpha,alpha)^2+
         dww.fun(w,q)*wgg.fun(y,x.gamma,gamma,x.alpha,alpha)
  return(result)
}
dwga.fun=function(w,q,y,x.gamma,gamma,x.alpha,alpha){
  result=dwww.fun(w,q)*wg.fun(y,x.gamma,gamma,x.alpha,alpha)*wa.fun(y,x.gamma,gamma,x.alpha,alpha)+
         dww.fun(w,q)*wga.fun(y,x.gamma,gamma,x.alpha,alpha)
  return(result)
}
dwgq.fun=function(w,q,y,x.gamma,gamma,x.alpha,alpha){
  result=dwwq.fun(w,q)*wg.fun(y,x.gamma,gamma,x.alpha,alpha)
  return(result)
}
dwaa.fun=function(w,q,y,x.gamma,gamma,x.alpha,alpha){
  result=dwww.fun(w,q)*wa.fun(y,x.gamma,gamma,x.alpha,alpha)^2+
         dww.fun(w,q)*waa.fun(y,x.gamma,gamma,x.alpha,alpha) 
  return(result)
}
dwaq.fun=function(w,q,y,x.gamma,gamma,x.alpha,alpha){
  result=dwwq.fun(w,q)*wa.fun(y,x.gamma,gamma,x.alpha,alpha)
  return(result)
}
dwqq.fun=function(w,q){
  result=dglgdqq(w,q)
  return(result)
}
###
Swww.fun=function(w,q){
  result=-dww.fun(w,q)
  return(result)
}
Swwq.fun=function(w,q){
  result=sglgdwq(w,q)
  return(result)
}
Swgg.fun=function(w,q,y,x.gamma,gamma,x.alpha,alpha){
  result=Swww.fun(w,q)*wg.fun(y,x.gamma,gamma,x.alpha,alpha)^2+
         Sww.fun(w,q)*wgg.fun(y,x.gamma,gamma,x.alpha,alpha)
  return(result)
}
Swga.fun=function(w,q,y,x.gamma,gamma,x.alpha,alpha){
  result=Swww.fun(w,q)*wg.fun(y,x.gamma,gamma,x.alpha,alpha)*wa.fun(y,x.gamma,gamma,x.alpha,alpha)+
         Sww.fun(w,q)*wga.fun(y,x.gamma,gamma,x.alpha,alpha)
  return(result)
}
Swgq.fun=function(w,q,y,x.gamma,gamma,x.alpha,alpha){
  result=Swwq.fun(w,q)*wg.fun(y,x.gamma,gamma,x.alpha,alpha)
  return(result)
}
Swaa.fun=function(w,q,y,x.gamma,gamma,x.alpha,alpha){
  result=Swww.fun(w,q)*wa.fun(y,x.gamma,gamma,x.alpha,alpha)^2+
         Sww.fun(w,q)*waa.fun(y,x.gamma,gamma,x.alpha,alpha)
  return(result)
}
Swaq.fun=function(w,q,y,x.gamma,gamma,x.alpha,alpha){
  result=Swwq.fun(w,q)*wa.fun(y,x.gamma,gamma,x.alpha,alpha)
  return(result)
}
Swqq.fun=function(w,q){
  result=sglgdqq(w,q)
  return(result)
}

######################
### Other function ###
######################
moments.GGD=function(gamma,alpha,q,org.time=0){
  sigma=exp(alpha)
  median.GGD=org.time+exp(gamma)*exp(sigma*qglgd(0.5,q))
  mean.GGD=org.time+exp(gamma)*(q^2)^(sigma/q)*gamma(q^(-2)+sigma/q)/gamma(q^(-2))
  var.GGD=mean.GGD^2*(gamma(q^(-2)+2*sigma/q)*gamma(q^(-2))/gamma(q^(-2)+sigma/q)^2-1)
  std.GGD=sqrt(var.GGD)

  return(list(median.GGD=median.GGD,mean.GGD=mean.GGD,std.GGD=std.GGD))
}
### Example
if(F){
  moments.GGD(3.751,-1.207,4.536,20)
}
