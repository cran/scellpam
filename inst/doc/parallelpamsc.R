## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(scellpam)

## -----------------------------------------------------------------------------
# Initially, state of debug is FALSE. Turn it on exclusively for the
# parallelpam part with
ScellpamSetDebug(FALSE,debparpam=TRUE)
# There is another parameter, debjmat, to turn on messages about
# binary matrix creation/manipulation. By default is FALSE but turn it on
# if you like with
# ScellpamSetDebug(FALSE,debparpam=TRUE,debjmat=TRUE)

## -----------------------------------------------------------------------------
# Create the matrix with row names V1 to V5000 and column names d1 to d500
nvec<-5000
ndim<-500
P<-matrix(runif(nvec*ndim),nrow=nvec)
rownames(P)<-paste0("V",1:nvec)
colnames(P)<-paste0("d",1:ndim)
# Write it to disk as a binary file in jmatrix format. Please,
# see vignette jmatrixsc.
JWriteBin(P,"datatest.bin",dtype="float",dmtype="full",
          comment="Synthetic problem data to test PAM")

## -----------------------------------------------------------------------------
JMatInfo("datatest.bin")

## -----------------------------------------------------------------------------
CalcAndWriteDissimilarityMatrix("datatest.bin","datatestL2.bin",
                                distype="L2",restype="float",
                                comment="L2 distance for vectors in
 jmatrix file datatest.bin",nthreads=0)

## -----------------------------------------------------------------------------
JMatInfo("datatestL2.bin")

## -----------------------------------------------------------------------------
L=ApplyPAM("datatestL2.bin",k=5,init_method="BUILD",max_iter=1000,
           nthreads=0)

## -----------------------------------------------------------------------------
Lbuild=ApplyPAM("datatestL2.bin",k=5,init_method="BUILD",max_iter=0)
Llab1=ApplyPAM("datatestL2.bin",k=5,init_method="LAB",max_iter=0)
Llab2=ApplyPAM("datatestL2.bin",k=5,init_method="LAB",max_iter=0)

## -----------------------------------------------------------------------------
Llab2Final=ApplyPAM("datatestL2.bin",k=5,init_method="PREV",
                    initial_med=Llab2$med)

## -----------------------------------------------------------------------------
# Which are the indexes of the points chosen as medoids?
L$med
#
# In which class has point 147 been classified?
L$clasif[147]
#
# And which is the index (row in the dissimilarity matrix)
# of the medoid closest to point 147?
L$med[L$clasif[147]]

## -----------------------------------------------------------------------------
min(L$clasif)
max(L$clasif)

## -----------------------------------------------------------------------------
S=CalculateSilhouette(Llab2$clasif,"datatestL2.bin",nthreads=0)

## -----------------------------------------------------------------------------
Sclus <- NumSilToClusterSil(Llab2$clasif,S)
library(cluster)
plot(Sclus)

## -----------------------------------------------------------------------------
Lfilt=FilterBySilhouetteQuantile(S,Llab2,"datatest.bin",
                                 "datatestFilt.bin","datatestL2.bin",
                                 "datatestL2Filt.bin",0.2)

## -----------------------------------------------------------------------------
Lfinal=ApplyPAM("datatestL2Filt.bin",k=length(Lfilt$med),
                init_method="PREV",initial_med=Lfilt$med)

## ---- results='hide'----------------------------------------------------------
d = GetSubdiag("datatestL2.bin")

## ---- eval=FALSE--------------------------------------------------------------
#  library(cluster)
#  clusterpam = pam(d,diss=TRUE,k=5)
#  print(sort(clusterpam$id.med))
#  print(sort(L$med))

## ---- eval=FALSE--------------------------------------------------------------
#  # Be patient, this may take some time...
#  Dm = GetJManyRows("datatestL2.bin",seq(1:nvec))

## ---- eval=FALSE--------------------------------------------------------------
#  library(ClusterR)
#  ClusterRpam = Cluster_Medoids(Dm,clusters=5)
#  print(sort(ClusterRpam$medoid_indices))
#  print(sort(L$med))

## ---- eval=FALSE--------------------------------------------------------------
#  TDparallelpam = GetTD(L,"datatestL2.bin")
#  
#  # This is to adapt cluster package output format to ours, since this is what our GetTD function expects...
#  Lcl = list()
#  Lcl$med = clusterpam$id.med
#  Lcl$clasif = clusterpam$clustering
#  TDcluster = GetTD(Lcl,"datatestL2.bin")
#  
#  # The same with ClusterR package:
#  LclR = list()
#  LclR$med = ClusterRpam$medoid_indices
#  LclR$clasif = ClusterRpam$clusters
#  TDClusterR = GetTD(LclR,"datatestL2.bin")

