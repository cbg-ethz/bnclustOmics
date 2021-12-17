library(bnclustOmics) #github
library(BiDAG) #on CRAN
#first the data object is needed
#a list of matrices, 1 per omics type
#rows are samples, columns are genes (features)
dataHCCt<-readRDS("testdata.rds")

#five omics
names(dataHCCt)

dim(dataHCCt$M)
dim(dataHCCt$PP)
#...etc

head(dataHCCt$M)
head(dataHCCt$PP)
#...etc

#mappings are needed if the IDs are different for different omics types
#so for example ENSEMBLE is often use for transcriptome and UNIPROT for proteome
#just one columnt is needed 'gene', the rownames have to be similar to IDs
#used for corresponding omics type

#create mappings object
mappings<-list()
mappings[["M"]]<-readRDS(("mappings/mutMAP.rds"))
mappings[["CN"]]<-readRDS(("mappings/CNAMAP.rds"))
mappings[["T"]]<-readRDS(("mappings/transcriptMAP.rds"))
mappings[["P"]]<-readRDS(("mappings/protMAP.rds"))
mappings[["PP"]]<-readRDS(("mappings/pprotMAP.rds"))

#examples, all columns apart from 'gene' are disregarded
head(mappings[["M"]])
head(mappings[["P"]])

#create omics object containing all neccessary mappings and other information needed for bnclustOmics
bnnames<-bnInfo(dataHCCt,c("b","o","c","c","c"),c("M","CN","T","P","PP"),
                 mappings)

bnnames$nb #number or M features
bnnames$totn #total number of features, all omics types
bnnames$allnamesonebn #names of each feature
#etc...

#initialize blacklist (optional)
blt<-blInit(bnnames,intra=c("T"),interXX=list(from=c("P"),to=c("T")),
            interXY=list(from=c("CN","CN","CN"),to=c("T","PP","P")))


#read the file containing interactions (prior information)
stringint<-readRDS("stringint.rds")
head(stringint)

#initialize penalization matrix (optional)
pmt<-bnclustOmics::penInit(bnnames,pfbase=2,intpf=1,intlist=stringint,intsame = 1, usescore=FALSE)
#2-2*interactions_score,
hist(pmt)
#blt<-pmt-1 #blacklist all non-string interactions

#run bnclustOmics
#if the number of nodes >100, better to be run on Euler
#for example, clustering for kclust=3 and total number of features ~770 it takes ca.20 hours
bnres<-bnclustOmics::bnclustOmics(dataHCCt,bnnames,
                                  blacklist=blt,
                                  edgepmat = pmt, seed=100, baseprob=0.6,
                                  maxEM=4, kclust=2, epmatrix = FALSE,
                                  startpoint = "mclustPCA",
                                  hardlim=6, deltahl=2,
                                  commonspace=TRUE)
#cluster membership
bnres$memb
#consensus membership, most likely the same
bnres$consmemb
bnres$consmemb==bnres$memb


#run for kclust=1,2,3,4,5...compare AICs, BICs, pick lowest
bnres$AIC #lower sample size
bnres$BIC #bigger sample size

#look at MAP networks
length(bnres$DAGs)
head(bnres$DAGs[[1]])
compareDAGs(bnres$DAGs[[1]],bnres$DAGs[[2]])

#look at consensus networks
cons01<-getModels(bnres,0.1) #getting consensus graphs with a posterior threshold of 0.1
cons05<-getModels(bnres,0.5) #getting consensus graphs with a posterior threshold of 0.5
cons09<-getModels(bnres,0.9) #getting consensus graphs with a posterior threshold of 0.9

compareDAGs(cons01[[1]],cons01[[2]])
compareDAGs(cons05[[1]],cons05[[2]])

#annotate all edges
allInteractions<-annotateEdges(bnres,bnnames,sump=1.2,minp=0.5,minkp=0.9,dblist=stringint)
head(allInteractions)
nrow(allInteractions)
length(which(allInteractions$db))
length(which(allInteractions$gene1==allInteractions$gene2))
#add correlations: columns cl1, cl2, etc...
allInteractions<-addCors(dataHCCt,allInteractions,bnres$memb)
head(allInteractions)

#21 out of 86 discovered interactions can be found in the string database, pen factor 2
#25 out of 77 interactions can be found in the string database, pen factor 4
#13 out of 29 interactions can be found in the string database, pen factor 100
#here also check for interactions between same genes, they are not part of DBs usually






