## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(scellpam)

## -----------------------------------------------------------------------------
# Initially, state of debug is FALSE. Turn it on with
ScellpamSetDebug(TRUE,debparpam=TRUE,debjmat=TRUE)
# We have also turned on debugging of the parallel PAM algorithm and of the binary matrix creation/manipulation
# but if this annoys you their default value is FALSE for both cases even when the scellpam debug is set to true, so just do
# ScellpamSetDebug(TRUE)

## ---- eval=FALSE--------------------------------------------------------------
#  CsvToJMat("countsfile.csv","countsfile.bin")

## ---- eval=FALSE--------------------------------------------------------------
#  CsvToJMat("countsfile.csv","countsfile.bin",mtype="sparse",csep=",",
#            ctype="raw",valuetype="float",transpose=FALSE,comment="")

## ---- eval=FALSE--------------------------------------------------------------
#  CsvToJMat("countsfile.csv","countsfile.bin",ctype="log1n",transpose=TRUE,
#              comment="Obtained from countsfile.csv")

## ---- eval=FALSE--------------------------------------------------------------
#  JMatToCsv("countsfile.bin","countsback.csv",csep=",",withquotes=FALSE)

## ---- eval=FALSE--------------------------------------------------------------
#  q <- readRDS("yourdata.rds")
#  dgCMatToJMat(q@assays$RNA@counts,"countsfile.bin",transpose=TRUE,
#               comment="Obtained from yourdata.rds")

## ---- eval=FALSE--------------------------------------------------------------
#  dgCMatToJMat(q@raw.data,"countsfile.bin",transpose=TRUE,
#               comment="Obtained from other data")

## ---- eval=FALSE--------------------------------------------------------------
#  Gr=GetSeuratGroups(q)

## ---- eval=FALSE--------------------------------------------------------------
#  suppressPackageStartupMessages({
#     library(splatter)
#     library(scater)
#  })
#  set.seed(1)
#  Splattersce <- mockSCE()
#  SceToJMat(Splattersce@assays@data@listData$counts,"mockSCE.bin",
#             mtype="full",ctype="log1n",transpose=TRUE,
#             comment="Generated by splatter with mockSCE() and normalized
#                      to log2(counts+1).")

## ---- eval=FALSE--------------------------------------------------------------
#  suppressPackageStartupMessages({
#     library(DuoClustering2018)
#  })
#  KuTCCsce <- sce_filteredExpr10_KumarTCC()
#  SceToJMat(KuTCCsce@assays$data@listData$counts,"KuTCC.bin",
#            mtype="full",ctype="log1n",transpose=TRUE,
#            comment="Generated by DuoCLustering with Kumar data and
#                  normalized to log2(counts+1).")

## -----------------------------------------------------------------------------
tmpfile=paste0(tempdir(),"/Trapnell.bin")
CsvToJMat("extdata/conquer_GSE52529_Trapnell_sample.csv",tmpfile,
          transpose=TRUE,comment="Experiment conquer GSE52529-GPL16791")
JMatInfo(tmpfile)

## -----------------------------------------------------------------------------
tmpfile=paste0(tempdir(),"/Trapnell.bin")
tmptxtfile=paste0(tempdir(),"/TrapnellInfo.txt")
JMatInfo(tmpfile,tmptxtfile)

## -----------------------------------------------------------------------------
v<-as.vector(read.table("extdata/Trapnell_remain_genes.csv")[[1]])
tmpfile1=paste0(tempdir(),"/Trapnell.bin")
tmpfiltfile1=paste0(tempdir(),"/Trapnell_filtered.bin")
FilterJMatByName(tmpfile1,v,tmpfiltfile1,namesat="cols")

## ---- eval=FALSE--------------------------------------------------------------
#  tmpfile1=paste0(tempdir(),"/Trapnell.bin")
#  tmppeafile1=paste0(tempdir(),"/TrapnellPearson.bin")
#  CalcAndWriteDissimilarityMatrix(tmpfile1,tmppeafile1,distype="Pearson",nthreads=0,
#                  comment="Pearson dissimilarities for cell in experiment conquer GSE52529-GPL16791")

## -----------------------------------------------------------------------------
tmpfile1=paste0(tempdir(),"/Trapnell.bin")
tmppeafile1=paste0(tempdir(),"/TrapnellPearson.bin")
CalcAndWriteDissimilarityMatrix(tmpfile1,tmppeafile1,nthreads=-1,
                                comment="Pearson dissimilarities for cell
                                         in experiment conquer
                                         GSE52529-GPL16791")

## -----------------------------------------------------------------------------
JMatInfo(tmppeafile1)

## -----------------------------------------------------------------------------
tmppeafile1=paste0(tempdir(),"/TrapnellPearson.bin")
L=ApplyPAM(tmppeafile1,k=10,init_method="BUILD",
           max_iter=1000,nthreads=-1)

## -----------------------------------------------------------------------------
tmppeafile1=paste0(tempdir(),"/TrapnellPearson.bin")
Lbuild=ApplyPAM(tmppeafile1,k=10,init_method="BUILD",max_iter=0)
Llab1=ApplyPAM(tmppeafile1,k=10,init_method="LAB",max_iter=0)
Llab2=ApplyPAM(tmppeafile1,k=10,init_method="LAB",max_iter=0)

## -----------------------------------------------------------------------------
tmppeafile1=paste0(tempdir(),"/TrapnellPearson.bin")
Llab2Final=ApplyPAM(tmppeafile1,k=10,init_method="PREV",
                    initial_med=Llab2$med,nthreads=-1)

## -----------------------------------------------------------------------------
# Which are the indexes of the points chosen as medoids?
L$med
#
# In which class has point 147 been classified?
L$clasif[147]
#
# And which is the index (row in the dissimilarity matrix) of the medoid closest to point 147?
L$med[L$clasif[147]]

## -----------------------------------------------------------------------------
min(L$clasif)
max(L$clasif)

## -----------------------------------------------------------------------------
tmppeafile1=paste0(tempdir(),"/TrapnellPearson.bin")
F <- ClassifAsDataFrame(L,tmppeafile1)

## -----------------------------------------------------------------------------
# Extract column 1 (cell name) of those rows for which distance to
# closest medoid (column 3) is 0
F[which(F[,3]==0),1]

## ---- eval=FALSE--------------------------------------------------------------
#  M=BuildAbundanceMatrix(L$clasif,Gr)

## -----------------------------------------------------------------------------
tmppeafile1=paste0(tempdir(),"/TrapnellPearson.bin")
S=CalculateSilhouette(Llab2$clasif,tmppeafile1,nthreads=-1)

## -----------------------------------------------------------------------------
Sclus <- NumSilToClusterSil(Llab2$clasif,S)
library(cluster)
plot(Sclus)

## -----------------------------------------------------------------------------
tmpfile1=paste0(tempdir(),"/Trapnell.bin")
tmpfiltfile1=paste0(tempdir(),"/TrapnellFilt.bin")
tmppeafile1=paste0(tempdir(),"/TrapnellPearson.bin")
tmppeafiltfile1=paste0(tempdir(),"/TrapnellPearsonFilt.bin")
Lfilt=FilterBySilhouetteQuantile(S,Llab2,tmpfile1,tmpfiltfile1,
                                 tmppeafile1,tmppeafiltfile1,0.2)

## -----------------------------------------------------------------------------
tmppeafiltfile1=paste0(tempdir(),"/TrapnellPearsonFilt.bin")
Ffilt=ClassifAsDataFrame(Lfilt,tmppeafiltfile1)

## -----------------------------------------------------------------------------
tmpfiltfile1=paste0(tempdir(),"/TrapnellFilt.bin")
tmpfiltcsvfile1=paste0(tempdir(),"/TrapnellFilt.csv")
JMatToCsv(tmpfiltfile1,tmpfiltcsvfile1)

## -----------------------------------------------------------------------------
tmppeafiltfile1=paste0(tempdir(),"/TrapnellPearsonFilt.bin")
Lfinal=ApplyPAM(tmppeafiltfile1,k=length(Lfilt$med),init_method="PREV",initial_med=Lfilt$med,nthreads=-1)

