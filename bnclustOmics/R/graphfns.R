#'Deriving consensus graphs
#'
#'This function derives a consensus grpah for each cluster according to a specified posterior probability.
#'@param bnres object of class 'bnclustOmics'
#'@param p posterior probability threshold
#'@return a list of adjacency matrices corresponding to consensus graphs representing discovered clusters
#'@author Polina Suter
#'@export
getModels<-function(bnres,p) {
  return(lapply(bnres$ep,getmodel,p))
}

getmodel<-function(ep, p) {
  n<-ncol(ep)
  labels<-colnames(ep)
  incidence <- matrix(rep(0, n * n), nrow = n, ncol = n)
  colnames(incidence)<-rownames(incidence)<-labels
  incidence[which(ep > p)] <- 1
  return(incidence)
}
