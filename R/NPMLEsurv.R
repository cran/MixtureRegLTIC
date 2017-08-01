#########################
### function: NPMLEsurv()
#########################
NPMLEsurv=function(formula,var.entry,var.weight=NULL,data,time.origin=0){

  cat("Running ... \n")

  dataset=data

  vars=split.vars(formula)
  response=vars$y

  covariates=list(names=sort(unique(c(vars$main,vars$i2,vars$i3))))

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
    index=which(dataset[,var.centime[3]]==1)
    if(length(index)>0){
      dataset[index,var.centime[2]]=dataset[index,var.centime[1]]
    }
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
    index=which(is.na(dataset[,var.entry]))
    if(length(index)>0) dataset[index,"status.entry"]=0
  }else{
    var.truntime=NULL
  }

  tolerance=1e-6

  est=turnbullFrydman(data=dataset,var.centime=var.centime,var.truntime=var.truntime,
                  covariates=covariates,time.origin=time.origin,var.weight=var.weight,tolerance=tolerance)

  est$method="NPMLE"

  cat(est$run.time); cat("\n")

  class(est)="NPMLEsurv"

  return(est)
}
