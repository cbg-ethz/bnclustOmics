#'Annotating edges from discovered networks
#'
#'This function makes a data frame which contains all pairs of nodes connected in cluster-specific networks
#'@param bnres an object of class 'bnclustOmics'; see \link{bnclustOmics}
#'@param bnnames an object of class 'bnInfo'; see \link{bnInfo}
#'@param sump threshold for the sum of posterior probabilities in all discovered networks
#'@param minp threshold for the minimum posterior probability in at least one network, when the sum of posterioirs is bigger than sump
#'@param minkp threshold for the minimum posterior probability in at least one network, when the sum of posterioirs is less than sump
#'@param maxkp (optional) threshold for the maximum posterior probability in at least one network; used to esclude cluster specific edges from the edges with high sum of posterioirs (>sump)
#'@param dblist a list of known interactions, discovered edges will be annotated is the edge is present in this list; two columns must be present 'gene1' and 'gene2'
#'@export
annotateEdges<-function(bnres,bnnames,sump=1.2,minp=0.5,minkp=0.9,
                        maxkp=NULL,dblist=NULL) {
  consnets<-getModels(bnres,min(minp,minkp))

  consp<-lapply(bnres$ep, getmodel, p=minp)
  comcons<-Reduce('+',bnres$ep)
  constot<-1*Reduce('|',consp)
  constot[which(comcons<sump)]<-0
  if(!is.null(maxkp)) {
    consmaxkp<-lapply(bnres$ep, getmodel, p=maxkp)
    constotmaxkp<-1*Reduce('|',consmaxkp)
    constot[which(constotmaxkp>0)]<-0
  }

  conskp<-lapply(bnres$ep, getmodel, p=minkp)
  constotkp<-1*Reduce('|',conskp)
  constotkp[which(comcons>sump)]<-0

  constot<-1*(constot|constotkp)
  intconstot<-getAllInt(bnnames,constot)
  intconstot<-checkInt(intconstot,dblist,"db")
  intconstot<-addpost(intconstot,bnres)


return(intconstot)
}
getAllInt<-function(bnnames,dag,bnparam=NULL){
  allnames<-rownames(dag)
  alltypes<-rep(bnnames$omics,bnnames$ns)
  names(alltypes)<-bnnames$allnamesonebn
  int<-NULL
  for(i in 1:nrow(dag)) {
    chs<-which(dag[i,]>0)
    if(length(chs)>0) {
      par<-allnames[i]
      vars<-names(chs)
      chs<-as.vector(chs)
      typepar<-as.vector(alltypes[par])
      typechs<-as.vector(alltypes[vars])
      int<-rbind(int,data.frame(from=par,to=vars,type1=typepar,type2=typechs))
    }
  }

  int<-transformInt(int,bnnames)
  rownames(int)<-c(1:nrow(int))
  return(int)
}
#'Adding cluster-specific correlations to annotated list of interactions
#'
#'This function adds cluster-specific correlations to the list of annotated interactions; see \link{annotateEdges}
#'@param omicdata a list of matrices orresponding to omics types: "M", "CN", "T", "P" and "PP"; at least one continuous type must be present
#'@param intlist annotated list of interactions produced by the function \link{annotateEdges}
#'@param memb vector of memberships learned by the function \link{bnclustOmics}
#'@export
addCors<-function(omicdata,intlist,memb) {
  addco<-NULL
  kclust<-length(unique(memb))
  for(i in 1:nrow(intlist)) {
    if(grepl(".CN",intlist$from[i])) newfrom<-sub("\\.CN","",intlist$from[i]) else newfrom<-intlist$from[i]
    addco_local<-c()
    for(j in 1:kclust) {
      addco_local<-c(addco_local,getcors_local(omicdata,newfrom,intlist$to[i],intlist$type1[i],intlist$type2[i],memb,j))
    }
    addco<-rbind(addco,addco_local)
  }
  colnames(addco)<-paste("cl",1:kclust,sep="")
  return(cbind(intlist,addco))
}
paramInt<-function(int,bnparam) {
  addval<-NULL
  for(i in 1:nrow(int)) {
    ch<-int$to[i]
    sumlocal<-summary(bnparam[[ch]])
    if(int$from[i]%in%rownames(sumlocal$coefficients)) {
      p<-sumlocal$coefficients[int$from[i],4]
      r2<-sumlocal$r.squared
      addval<-rbind(addval, data.frame(pval=p,r2=r2))
    } else {
      addval<-rbind(addval, data.frame(pval=1,r2=0))
    }
  }
  return(cbind(int,addval))
}
