/*
 *
 * Copyright (C) 2023 Juan Domingo (Juan.Domingo@uv.es)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <dissimmat.h>

extern unsigned char DEB;

//' CalcAndWriteDissimilarityMatrix
//'
//' Writes a binary symmetric matrix with the dissimilarities between ROWS of the data stored in a binary matrix in the scellpam package format.\cr
//' Notice that, differently from the common practice in single cell, the rows represent cells. This is for efficiency reasons and it is transparent
//' to the user, as long as he/she has generated the binary matrix (with CsvToBinMat, dgCMatToBinMat or SceToBinMat) using the option transpose=TRUE.\cr
//' The input matrix of vectors can be a full or a sparse matrix. Output matrix type can be float or double type (but look at the comments in 'Details').
//'
//' The parameter restype forces the output to be a matrix of either floats or doubles. Precision of float is normally good enough; but if you need
//' double precision (may be because you expect your results to be in a large range, two to three orders of magnitude), change it.\cr
//' Nevertheless, notice that this at the expense of double memory usage, which is QUADRATIC with the number of individuals (rows) in your input matrix.
//'
//' @param ifname   A string with the name of the file containing the counts as a binary matrix, as written by CsvToBinMat, dgCMatToBinMat or SceToBinMat
//' @param ofname   A string with the name of the binary output file to contain the symmetric dissimilarity matrix.
//' @param distype  The dissimilarity to be calculated. It must be one of these strings: 'L1', 'L2', 'Pearson', 'Cos' or 'WEuc'.\cr
//'                 Respectively: L1 (Manhattan), L2 (Euclidean), Pearson (Pearson dissimilarity), Cos (cosine distance), WEuc (weigthed Euclidean, with inverse-stdevs as weights).\cr
//'                 Default: 'L2'.
//' @param restype  The data type of the result. It can be one of the strings 'float' or 'double'. Default: float (and don't change it unless you REALLY need to...).
//' @param comment  Comment to be added to the dissimilary matrix. Default: "" (no comment)
//' @param nthreads Number of threads to be used for the parallel calculations with this meaning:\cr
//'                 -1: don't use threads.\cr
//'                  0: let the function choose according to the number of individuals (cells) and to the number of available cores.\cr
//'                  Any possitive number > 1: use that number of threads. You can use even more than cores, but this is discouraged and raises a warning.\cr
//'                 Default: 0.
//' @return   No return value, called for side effects (creates a file)
//' @examples
//' Rf <- matrix(runif(50000),nrow=100)
//' tmpfile1=paste0(tempdir(),"/Rfullfloat.bin")
//' JWriteBin(Rf,tmpfile1,dtype="float",dmtype="full",
//'           comment="Full matrix of floats, 100 rows, 500 columns")
//' JMatInfo(tmpfile1)
//' tmpdisfile1=paste0(tempdir(),"/RfullfloatDis.bin")
//' # Distance file calculated from the matrix stored as full
//' CalcAndWriteDissimilarityMatrix(tmpfile1,tmpdisfile1,distype="L2",
//'                          restype="float",comment="L2 distance matrix from full",nthreads=0)
//' JMatInfo(tmpdisfile1)
//' tmpfile2=paste0(tempdir(),"/Rsparsefloat.bin")
//' JWriteBin(Rf,tmpfile2,dtype="float",dmtype="sparse",
//'                          comment="Sparse matrix of floats, 100 rows, 500 columns")
//' JMatInfo(tmpfile2)
//' # Distance file calculated from the matrix stored as sparse
//' tmpdisfile2=paste0(tempdir(),"/RsparsefloatDis.bin")
//' CalcAndWriteDissimilarityMatrix(tmpfile2,tmpdisfile2,distype="L2",
//'                          restype="float",comment="L2 distance matrix from sparse",nthreads=0)
//' JMatInfo(tmpdisfile2)
//' # Read both versions
//' Dfu<-GetJManyRows(tmpdisfile1,c(1:nrow(Rf)))
//' Dsp<-GetJManyRows(tmpdisfile2,c(1:nrow(Rf)))
//' # and compare them
//' max(Dfu-Dsp)
//' @export
// [[Rcpp::export]]
void CalcAndWriteDissimilarityMatrix(std::string ifname, std::string ofname, std::string distype="L2", std::string restype="float", std::string comment="",int nthreads=0)
{
 if ((distype != "L1") && (distype != "L2") && (distype != "Pearson") && (distype != "Cos") && (distype != "WEuc"))
 {
  Rcpp::stop("Parameter distype must be one of 'L1', 'L2', 'Pearson', 'Cos' or 'WEuc'.\n");
  return;
 }
 if ((restype != "float") && (restype != "double"))
 {
  Rcpp::stop("Parameter restype must be one of 'float' or 'double'.\n");
 }
 
 unsigned char dtype;
 if (distype=="L1")
  dtype=DL1;
 if (distype=="L2")
  dtype=DL2;
 if (distype=="Pearson")
  dtype=DPe;
 if (distype=="Cos")
  dtype=DCo;
 if (distype=="WEuc")
  dtype=DWe; 
  
 unsigned char mt,ct,e,md;
 indextype nr,nc;
 MatrixType(ifname,mt,ct,e,md,nr,nc);
 
 if (DEB & DEBPP)
  Rcpp::Rcout << "Input matrix is ";
 switch (mt)
 {
  case MTYPEFULL: if (DEB & DEBPP) 
                   Rcpp::Rcout << "a full matrix ";
                  break;
  case MTYPESPARSE: if (DEB & DEBPP)
                     Rcpp::Rcout << "a sparse matrix ";
                    break;
  case MTYPESYMMETRIC: if (DEB & DEBPP)
                        Rcpp::Rcout << "a symmetric matrix. This is not allowed; it must be full or sparse.\n";
                       Rcpp::stop("Invalid matrix type.\n");
                       break;
  default: if (DEB & DEBPP)
            Rcpp::Rcout << "of unknown type (neither full, sparse of symmetric). Was it created with jmatrix?\n";
           Rcpp::stop("Unknown matrix type.\n");
           break;
 }
 switch (ct)
 {
  case FTYPE: if (DEB & DEBPP)
               Rcpp::Rcout << " with elements of type 'float' and size (" << nr << "," << nc << ")\n";
              break;
  case DTYPE: if (DEB & DEBPP)
               Rcpp::Rcout << " with elements of type 'double' and size (" << nr << "," << nc << ")\n";
              break;
  // @TODO: other types of input matrix might be allowed in successive versions. To de done.
  default: if (DEB & DEBPP)
            Rcpp::Rcout << " with elements which are neither 'float' nor 'double'. This is not allowed to calculate dissimilarity matrix. Sorry.\n";
           Rcpp::stop("Data type of input matrix not allowed.\n");
           break;
 }
 
 unsigned int nt=ChooseNumThreads(nthreads);
 
 switch (mt)
 {
  case MTYPEFULL: 
  {
   if (restype=="float")
   {
    if (ct==FTYPE)
     CalcAndWriteAuxFull<float,float>(ifname,ofname,dtype,nt,comment);
    else
     CalcAndWriteAuxFull<double,float>(ifname,ofname,dtype,nt,comment);
   }
   else
   {
    if (ct==FTYPE)
     CalcAndWriteAuxFull<float,double>(ifname,ofname,dtype,nt,comment);
    else
     CalcAndWriteAuxFull<double,double>(ifname,ofname,dtype,nt,comment); 
   }
  }
  break;
  case MTYPESPARSE:
  {
   if (restype=="float")
   {
    if (ct==FTYPE)
     CalcAndWriteAuxSparse<float,float>(ifname,ofname,dtype,nt,comment);
    else
     CalcAndWriteAuxSparse<double,float>(ifname,ofname,dtype,nt,comment);
   }
   else
   {
    if (ct==FTYPE)
     CalcAndWriteAuxSparse<float,double>(ifname,ofname,dtype,nt,comment);
    else
     CalcAndWriteAuxSparse<double,double>(ifname,ofname,dtype,nt,comment);
   }
  }
  break;
  default: Rcpp::stop("Unknown error. Matrix type was supposed to have been checked before.\n"); break;
 }
 return;
}


