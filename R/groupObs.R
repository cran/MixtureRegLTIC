##########################
### function: labelsData()
##########################
labelsData=function(data,labels=list()){

  data=as.data.frame(data)
  nobs=nrow(data); nvar=ncol(data)
  for(i in 1:nvar){
    if(is.factor(data[,i])) data[,i]=as.character(data[,i])
  }

  for(i in nvar:1){
    data[]=data[order(data[,i],decreasing=F),]
  }

  if(length(labels)==nvar){
    num.match.labels=0
    for(i in 1:nvar){
      if(length(labels[[i]])==length(unique(data[[i]]))) num.match.labels=num.match.labels+1
    }
    if(num.match.labels==nvar){
      for(i in nvar:1){
        data.var.i=data[,i];
        levels.var.i=sort(unique(data.var.i),decreasing=F)
        for(j in 1:length(levels.var.i)){
          index=which(data.var.i==levels.var.i[j])
          data[index,i]=labels[[i]][j]
        }
      }
    }
  }

  labels.data=data[,1]
  if(nvar>=2){
    for(i in 2:nvar){
      labels.data=paste(labels.data,data[,i],sep=", ")
    }
  }

  return(unique(labels.data))
}
### Example: labelsData()
if(F){
  data=OrchardSprays[c(1,2,3,4)]
  labels=list(c("d2","d3","d4","d5","d6","d7","d8","d9","d10","d12","d13","d14","d15","d16","d17","d19","d20",
                "d22","d24","d27","d28","d29","d36","d39","d43","d44","d47","d51","d55","d57","d60","d61","d69",
                "d71","d72","d76","d77","d80","d81","d84","d86","d87","d90","d92","d95","d114","d127","d130"),
              c("r1","r2","r3","r4","r5","r6","r7","r8"),
              c("c1","c2","c3","c4","c5","c6","c7","c8"),
              c("tA","tB","tC","tD","tE","tF","tG","tH"))[c(1,2,3,4)]
  out=labelsData(data);print(out)
  out=labelsData(data,labels);print(out)
  out=labelsData(data[,4]);print(out)
}

########################
### function: groupObs()
########################
groupObs=function(data,labels=list()){

  data=as.data.frame(data)
  nobs=nrow(data); nvar=ncol(data)
  for(i in 1:nvar){
    if(is.factor(data[,i])) data[,i]=as.character(data[,i])
  }

  group=list()

  ### "group$labels"
  group$labels=labelsData(data)

  ### "group$levels"
  group$levels=seq(0,length(group$labels)-1)

  ### "group$levelsTrue"
  group$levelsTrue=rep(F,length(group$labels))

  ### "group$obs"
  group$obs=data[,1]
  if(nvar>=2){
    for(i in 2:nvar){
      group$obs=paste(group$obs,data[,i],sep=", ")
    }
  }
  group$obs=match(group$obs,group$labels)-1

  ### "group$levelsTrue"
  index=match(as.numeric(names(table(group$obs))),group$levels)
  if(length(index)>0) group$levelsTrue[index]=T

  ### "group$labels"
  if(length(labels)==nvar){
    group$labels=labelsData(data,labels)
  }
  group$labels=as.character(group$labels)

  return(group)
}
### Example: groupObs()
if(F){
  data=OrchardSprays[c(1,2,3,4)]
  labels=list(c("d2","d3","d4","d5","d6","d7","d8","d9","d10","d12","d13","d14","d15","d16","d17","d19","d20",
                "d22","d24","d27","d28","d29","d36","d39","d43","d44","d47","d51","d55","d57","d60","d61","d69",
                "d71","d72","d76","d77","d80","d81","d84","d86","d87","d90","d92","d95","d114","d127","d130"),
              c("r1","r2","r3","r4","r5","r6","r7","r8"),
              c("c1","c2","c3","c4","c5","c6","c7","c8"),
              c("tA","tB","tC","tD","tE","tF","tG","tH"))[c(1,2,3,4)]
  out=groupObs(data); print(out)
  out=groupObs(data,labels); print(out)
  out=groupObs(data[,3:4]); print(out)
}
