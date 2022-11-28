## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(scellpam)

## -----------------------------------------------------------------------------
# Initially, state of debug is FALSE. Turn it on exclusively for
# the jmatrix part of scellpam with
ScellpamSetDebug(FALSE,debparpam=FALSE,debjmat=TRUE)

## -----------------------------------------------------------------------------
# Create a 6x8 matrix of random values
Rf <- matrix(runif(48),nrow=6)
# Set row and column names for it
rownames(Rf) <- c("A","B","C","D","E","F")
colnames(Rf) <- c("a","b","c","d","e","f","g","h")
# Let's see the matrix
Rf
# and write it as the binary file Rfullfloat.bin
JWriteBin(Rf,"Rfullfloat.bin",dtype="float",dmtype="full",
          comment="Full matrix of floats")
# Also, you can write it with double data type:
JWriteBin(Rf,"Rfulldouble.bin",dtype="double",dmtype="full",
          comment="Full matrix of doubles")

## -----------------------------------------------------------------------------
# Information about the float binary file
JMatInfo("Rfullfloat.bin")
# Same information about the double binary file
JMatInfo("Rfulldouble.bin")

## -----------------------------------------------------------------------------
# Reads row 1 into vector vf. Float values inside the file are
# promoted to double.
(vf<-GetJRow("Rfullfloat.bin",1))

## -----------------------------------------------------------------------------
# Checks the precision lost
max(abs(Rf[1,]-vf))

## -----------------------------------------------------------------------------
vd<-GetJRow("Rfulldouble.bin",1)
max(abs(Rf[1,]-vd))

## -----------------------------------------------------------------------------
# Read column number 3
(vf<-GetJCol("Rfullfloat.bin",3))
# Test precision
max(abs(Rf[,3]-vf))
# Read row with name C
(vf<-GetJRowByName("Rfullfloat.bin","C"))
# Read column with name c
(vf<-GetJColByName("Rfullfloat.bin","c"))
# Get the names of all rows or columns as vectors of R strings
(rn<-GetJRowNames("Rfullfloat.bin"))
(cn<-GetJColNames("Rfullfloat.bin"))
# Get the names of rows and columns simultaneosuly as a list of two elements
(l<-GetJNames("Rfullfloat.bin"))
# Get several rows at once. The returned matrix has the rows in the
# same order as the passed list,
# and this list can contain even repeated values
(vm<-GetJManyRows("Rfullfloat.bin",c(1,4)))

# Of course, columns can be extrated equally
(vc<-GetJManyCols("Rfulldouble.bin",c(1,4)))
# and similar functions are provided for extracting by names:
(vm<-GetJManyRowsByNames("Rfulldouble.bin",c("A","D")))
(vc<-GetJManyColsByNames("Rfulldouble.bin",c("a","d")))

## -----------------------------------------------------------------------------
# Generation of a 6x8 sparse matrix
Rsp <- matrix(rep(0,48),nrow=6)
sparsity <- 0.1
nnz <- round(48*sparsity)
where <- floor(47*runif(nnz))
val <- runif(nnz)
for (i in 1:nnz)
{
 Rsp[floor(where[i]/8)+1,(where[i]%%8)+1] <- val[i]
}
rownames(Rsp) <- c("A","B","C","D","E","F")
colnames(Rsp) <- c("a","b","c","d","e","f","g","h")
# Let's see the matrix
Rsp
# Write the matrix as sparse with type float
JWriteBin(Rsp,"Rspafloat.bin",dtype="float",dmtype="sparse",
          comment="Sparse matrix of floats")

## -----------------------------------------------------------------------------
JMatInfo("Rspafloat.bin")

## -----------------------------------------------------------------------------
Rns <- matrix(runif(49),nrow=7)
Rsym <- 0.5*(Rns+t(Rns))
rownames(Rsym) <- c("A","B","C","D","E","F","G")
colnames(Rsym) <- c("a","b","c","d","e","f","g")
# Let's see the matrix
Rsym
# Write the matrix as symmetric with type float
JWriteBin(Rsym,"Rsymfloat.bin",dtype="float",dmtype="symmetric",
          comment="Symmetric matrix of floats")
# Get the information 
JMatInfo("Rsymfloat.bin")

## -----------------------------------------------------------------------------
(d<-GetSubdiag("Rsymfloat.bin"))

