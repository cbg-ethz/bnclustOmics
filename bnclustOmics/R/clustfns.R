# https://doi.org/10.1038/s41467-018-06867-x
# author: Kuipers et al.
propersample <- function(x){if(length(x)==1) x else sample(x,1)}
# https://doi.org/10.1038/s41467-018-06867-x
# author: Kuipers et al.
calcloglike <- function(samplescores,tau) {
  # samplescores<-t(samplescores)
  maxscorey<-apply(samplescores,1,max) # find the max of each column
  loglike<-sum(log(colSums(t(exp(samplescores-maxscorey))*tau))+maxscorey) # remove max for numerical stability and exponentiate
  return(loglike)
}
# https://doi.org/10.1038/s41467-018-06867-x
# author: Kuipers et al.
reassignsamples <- function(samplescores,numsamps){
  newclustermembership <-rep(0,numsamps) # to store the new cluster
  for(s in 1:numsamps){ # run through the samples
    clusterscores<-samplescores[s,]
    maxscorey<-max(clusterscores) # take the maximum
    maxscoreelem<-which(clusterscores==maxscorey)
    newclustermembership[s]<-propersample(maxscoreelem)
  }
  return(newclustermembership)
}
# https://doi.org/10.1038/s41467-018-06867-x
# author: Kuipers et al.
relativeprobs <- function(samplescores,numsamps){
  relativeprobabs <-rep(0,numsamps) # to store the relative probabilities
  for(s in 1:numsamps){ # run through the samples
    clusterscores<-samplescores[s,]
    maxscorey<-max(clusterscores) # take the maximum
    shifty<-exp(clusterscores-maxscorey)
    rescaley<-shifty/sum(shifty)
    relativeprobabs[s]<-max(rescaley) # relative probabilities
  }
  return(relativeprobabs)
}
# https://doi.org/10.1038/s41467-018-06867-x
# author: Kuipers et al.
allrelativeprobs <- function(samplescores,numsamps){
  relativeprobabs <-samplescores # to store the relative probabilities
  for(s in 1:numsamps){ # run through the samples
    clusterscores<-samplescores[s,]
    maxscorey<-max(clusterscores) # take the maximum
    shifty<-exp(clusterscores-maxscorey)
    rescaley<-shifty/sum(shifty)
    relativeprobabs[s,]<-rescaley # relative probabilities
  }
  return(relativeprobabs)
}
# https://doi.org/10.1038/s41467-018-06867-x
# author: Kuipers et al.
relativeprobswithtau <- function(sampleprobs,tau){
  temp<-tau*t(sampleprobs)
  relativeprobabswithtau<-1/colSums(temp)*t(temp) # ugly code
  return(relativeprobabswithtau)
}
# https://doi.org/10.1038/s41467-018-06867-x
# author: Kuipers et al.
avescore <- function(samplescores,numsamps){
  averagescores <-rep(0,numsamps) # to store the relative probabilities
  for(s in 1:numsamps){ # run through the samples
    clusterscores<-samplescores[s,]
    maxscorey<-max(clusterscores) # take the maximum
    shifty<-exp(clusterscores-maxscorey) # exponentiate
    rescaley<-log(mean(shifty))+maxscorey # mean and turn back to log
    averagescores[s]<-rescaley # averagescore
  }
  return(averagescores)
}
# https://doi.org/10.1038/s41467-018-06867-x
# author: Kuipers et al.
reassignsamplesprop <- function(samplescores,numsamps,gamma){
  newclustermembership <-rep(0,numsamps) # to store the new cluster
  for(s in 1:numsamps){ # run through the samples
    clusterscores<-samplescores[s,]*gamma
    maxscorey<-max(clusterscores) # take the maximum
    shifty<-exp(clusterscores-maxscorey)
    rescaley<-shifty/sum(shifty)
    scorelength<-length(rescaley)
    newclustermembership[s]<-sample.int(scorelength,1,prob=rescaley) # sample according to scores
  }
  return(newclustermembership)
}
# https://doi.org/10.1038/s41467-018-06867-x
# author: Kuipers et al.
comparelik<-function(assignprogress) {
  whichmax<-which.max(unlist(lapply(assignprogress,function(x)x$likel[length(x$likel)])))
}



#'Comparing estimated and ground truth membership
#'
#'This function compares ground truth cluster labels to estimated cluster labels.
#'
#'@param k number of clusters
#'@param truememb ground truth labels
#'@param estmemb estimated labels
#'@import stats
#'@import clue
#'@import BiDAG
#'@import mclust
#'@import RBGL
#'@import graph
#'@import gRbase
#'@export
checkmembership <- function(k,truememb,estmemb) {
  relabelmatrix<-matrix(nrow=k,ncol=k) #i row_num label, j col_num estmemb
  estlabels<-list()
  truelabels<-list()
  for (i in 1:k) {
    truelabels[[i]]<-which(truememb==i)
  }
  for (i in 1:k) {
    estlabels[[i]]<-which(estmemb==i)
  }
  for (i in 1:k) {
    for (j in 1:k) {
      relabelmatrix[i,j]<-length(which(estlabels[[i]]%in%truelabels[[j]]))
    }
  }
  rowcol<-solve_LSAP(relabelmatrix,maximum = TRUE)
  res<-list()
  res$relabel<-as.vector(rowcol)
  res$ncorr<-0
  for (j in 1:min(k,k)) {
    res$ncorr<-res$ncorr+relabelmatrix[j,res$relabel[j]]
  }
  res$relabelmatrix<-relabelmatrix
  return(res)
}
generatetriple<-function(n) {
  resmat<-matrix(nrow=n,ncol=3)
  for (i in 1:n) {
    ordr<-sample.int(3,3)
    res<-vector(length=3)
    res[1]<-runif(1,min=0,max=1)
    res[2]<-runif(1,min=0,max=1-res[1])
    res[3]<-1-res[1]-res[2]
    resmat[i,ordr]<-res
  }
  return(resmat)
}
generatefour<-function(n) {
  resmat<-matrix(nrow=n,ncol=4)
  for (i in 1:n) {
    ordr<-sample.int(4,4)
    res<-vector(length=4)
    res[1]<-runif(1,min=0,max=1)
    res[2]<-runif(1,min=0,max=1-res[1])
    res[3]<-runif(1,min=0,max=1-res[1]-res[2])
    res[4]<-1-res[1]-res[2]-res[3]
    resmat[i,ordr]<-res
  }
  return(resmat)
}
#changed baseprob
generatevec<-function(n,k,membvec=NULL,all1=FALSE) {
  resmat<-matrix(nrow=n,ncol=k)
  if(is.null(membvec)) {
    for (i in 1:n) {
      res<-sample(10:(15+k),k,replace=TRUE)
      res<-res/sum(res)
      resmat[i,]<-res
    }} else {
      if(!all1) {
        for(i in 1:n) {
          baseprob<-1/(k+2)
          probvec<-rep(baseprob,k)
          probvec[membvec[i]]<-3*baseprob
          resmat[i,]<-probvec
        }
      } else {
        for(i in 1:n) {
          probvec<-rep(0,k)
          probvec[membvec[i]]<-1
          resmat[i,]<-probvec
        }
      }
    }
  return(resmat)
}
generatevec2<-function(n,k,membvec=NULL,baseprob=3/(k+2)) {
  resmat<-matrix(nrow=n,ncol=k)
  if(is.null(membvec)) {
    for (i in 1:n) {
      res<-sample(10:(15+k),k,replace=TRUE)
      res<-res/sum(res)
      resmat[i,]<-res
    }} else {
      for(i in 1:n) {
        probvec<-rep((1-baseprob)/(k-1),k)
        probvec[membvec[i]]<-baseprob
        resmat[i,]<-probvec
      }
    }
  return(resmat)
}
newLabels<-function(relabvec) {
  changetovec<-vector()
  for(i in 1:(length(relabvec))) {
    changetovec[i]<-which(relabvec==i)
  }
  return(changetovec)
}
relabMembership<-function(estMemb, changeto) { #changeto 3,2,1,4
  newMemb<-vector()
  for(i in 1:length(changeto)) {
    whichy<-which(estMemb==i) #find all labels equalling 1
    newMemb[whichy]<-changeto[i] #change to 3, i.e. changeto[1]
  }
  return(newMemb)
}
relabLambdas<-function(estLambdas, changeto) {
  newLamb<-estLambdas
  for(i in 1:length(changeto)) {
    newLamb[,changeto[i]]<-estLambdas[,i]
  }
  return(newLamb)
}
relabDAGs<-function(estDAGs, changeto) {
  newDAGs<-estDAGs
  for(i in 1:length(changeto)) {
    newDAGs[[changeto[i]]]<-estDAGs[[i]]
  }
  return(newDAGs)
}



