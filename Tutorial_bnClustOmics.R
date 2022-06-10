#toydata, mappings, intlist

library(bnclustOmics) #CRAN
library(BiDAG) #on CRAN
#first the data object is needed
#a list of matrices, 1 per omics type
#rows are samples, columns are genes (features)
data(toydata)

dim(toydata$M)
head(toydata$M)

dim(toydata$P)
head(toydata$P)
#...etc

#mappings are needed if the IDs are different for different omics types
#so for example ENSEMBLE is often use for transcriptome and UNIPROT for proteome
#for co structing blacklists and penalization matrices it is important to pass mappings of all IDs to gene names
#just one column is needed 'gene', the rownames have to be similar to gene IDs withing each omic types
data(mappings)
head(mappings[["M"]])
head(mappings[["PP"]])


#create omics object containing all neccessary mappings and other information needed for bnclustOmics
bnnames<-bnInfo(toydata,c("b","o","c","c","c"),c("M","CN","T","P","PP"),
                 mappings)
bnnames

#initialize blacklist (optional)
#we blacklist the edges between variables of type "T" (transcriptome): intra=c("T")

#we blacklist the edges from variables representing gene X of type "P"
#to variables of type "T" representing the same gene: interXX=list(from=c("P"),to=c("T"))
#note that the edges in the other directions are allowed

#we blacklist the edges from variables representing gene X of type "CN"
#to variables of type "T", "P" and "PP" representing the gene Y:
#interXY=list(from=c("CN","CN","CN"),to=c("T","P","PP"))

blt<-blInit(bnnames,intra=c("T"),interXX=list(from=c("P"),to=c("T")),
            interXY=list(from=c("CN","CN","CN"),to=c("T","P","PP")))

#read the file containing interactions from the string database (prior information)
data(stringint)
head(stringint)

#initialize penalization matrix (optional)
#intpf=1, means we do not penalize interactions that are found in the database
#intsame=1, we do not penalize interactions between nodes representing the same gene but different omics types
#pfbase=2, we penalize edges by a factor of 2 if th are not found in the database
pmt<-bnclustOmics::penInit(bnnames,pfbase=2,intpf=1,intlist=stringint,intsame = 1, usescore=TRUE)
#2-2*interactions_score
hist(pmt)
#blt<-pmt-1 #blacklist all non-string interactions

#run bnclustOmics, 10 minutes for toy dataset, k=2
#if the number of nodes >100, it may take a while, >1hr
#for example, clustering for kclust=3 and total number of features ~770 it takes ca.20 hours
bnres<-bnclustOmics(toydata,bnnames,blacklist=blt,
                    edgepmat = pmt, kclust=2)

clusters(bnres)
all(clusters(bnres)==clusters(bnres,consensus=TRUE))

#run for kclust=1,2,3,4,5...compare AICs, BICs, pick lowest, see function chooseK()
bnres$AIC #better for lower sample size
bnres$BIC #better for larger sample size


#look at consensus networks
cons01<-getModels(bnres,0.1) #getting consensus graphs with a posterior threshold of 0.1
cons05<-getModels(bnres,0.5) #getting consensus graphs with a posterior threshold of 0.5

compareDAGs(cons01[[1]],cons01[[2]])
compareDAGs(cons05[[1]],cons05[[2]])

#annotate all edges
allInteractions<-annotateEdges(bnres,bnnames,sump=1.2,minp=0.5,minkp=0.9,dblist=stringint)
head(allInteractions)
nrow(allInteractions)
#number of interactions for in the database
length(which(allInteractions$db))
allInteractions[allInteractions$db,]
#number of interactions between nodes representing the same genes
length(which(allInteractions$gene1==allInteractions$gene2))
allInteractions[allInteractions$gene1==allInteractions$gene2,]

#plot node neighborhood
node<-"P15088"
plotNode(allInteractions,node)

#different threshold and font size
plotNode(localint,node,p=0.9,cex=0.7)


#different node
node<-"CTNNB1"
plotNode(allInteractions,node,p=0.9,cex=0.7)

#simulated data
bnnames<-bnInfo(simdata,c("b","c"),c("M","T"))
bnres<-list()
for(k in 2:4) {
bnres[[k]]<-bnclustOmics::bnclustOmics(simdata,bnnames,maxEM=4, kclust=k,
                                  startpoint = "mclustPCA")
}

#clustering accuracy
checkmembership(clusters(bnres[[2]]),simclusters)
checkmembership(clusters(bnres[[3]]),simclusters)
checkmembership(clusters(bnres[[4]]),simclusters)

#the optimal number of clusters is 3 according to both AIC and BIC
chooseK(bnres,fun="BIC")
chooseK(bnres,fun="AIC")

#to check the structure fit we first need to relabel according to
#corresponance between clustering labels
bnres[[3]]<-relabelSimulation(bnres[[3]],simclusters)

#compare MAP estimates to ground truth
compareDAGs(dags(bnres[[3]])[[1]],simdags[[1]])
compareDAGs(dags(bnres[[3]])[[2]],simdags[[2]])
compareDAGs(dags(bnres[[3]])[[3]],simdags[[3]])

#getting consensus models

#threshold of 0.5
cons05<-getModels(bnres[[3]],0.5)

#compare consensus estimates (p=0.5) to ground truth
#TPR is better, SHD is better than for MAP graphs
compareDAGs(cons05[[1]],simdags[[1]])
compareDAGs(cons05[[2]],simdags[[2]])
compareDAGs(cons05[[3]],simdags[[3]])

#threshold of 0.9
cons09<-getModels(bnres[[3]],0.9)

#compare consensus estimates (p=0.9) to ground truth
#TPR is worse, but SHD is better than for MAP graphs and consensus graphs with p=0.5
compareDAGs(cons09[[1]],simdags[[1]])
compareDAGs(cons09[[2]],simdags[[2]])
compareDAGs(cons09[[3]],simdags[[3]])

allInteractions<-annotateEdges(bnres[[3]],bnnames,sump=1.2,minp=0.5,minkp=0.9,dblist=simint)
node1<-names(table(allInteractions$from)[order(table(allInteractions$from),decreasing = TRUE)[1]])
node2<-names(table(allInteractions$to)[order(table(allInteractions$to),decreasing = TRUE)[1]])

#number of all annotated interactions
nrow(allInteractions)
#number of true positives
length(which(allInteractions$db==TRUE))

#plotting neighborhoods of node "T5"
plotNode(allInteractions,node1,p=0.5,cex=0.7,dbcheck=FALSE)
#check if interaction is in the ground truth graphs
plotNode(allInteractions,node1,p=0.5,cex=0.7,dbcheck=TRUE)


#plotting neighborhoods of node "T43"
plotNode(allInteractions,node2,p=0.5,cex=0.7,dbcheck=FALSE)
#check if interaction is in the ground truth graphs
plotNode(allInteractions,node2,p=0.5,cex=0.7,dbcheck=TRUE)






