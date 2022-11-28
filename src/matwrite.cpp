/*
 *
 * Copyright (C) 2022 Juan Domingo (Juan.Domingo@uv.es)
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
 
 
#include  <Rcpp.h>

#include <fullmatrix.h>
#include <sparsematrix.h>
#include <symmetricmatrix.h>
#include <cmath>

extern unsigned char DEB;

// ********** Functions to write a R matrix in binary format

// Auxiliary functions to write matrices in jmatrix binary format

template <typename T>
void WriteRMatrixAsBin(unsigned char mtype,std::string fname,Rcpp::NumericMatrix &M,std::string comment)
{
 indextype nr=indextype(M.nrow());
 indextype nc=indextype(M.ncol());
 if ((mtype==MTYPESYMMETRIC) && (nr!=nc))
  Rcpp::stop("Symmetric matrices must be square to be written in jmatrix binary format.\n");
  
 indextype nrl=0,ncl=0;
 
 Rcpp::StringVector rownames=Rcpp::StringVector();
 Rcpp::StringVector colnames=Rcpp::StringVector();
 
 if (M.hasAttribute("dimnames"))
 {
  Rcpp::List q4 = M.attr("dimnames");
  
  if (q4[0]!=R_NilValue)
  {
   rownames=q4[0];
   nrl=indextype(rownames.length());
  }
  else 
   nrl=0;
   
  if ((nrl>0) && (nrl!=nr))
   Rcpp::stop("Strange Matrix object. The number of rows in the matrix differs from the length of the vector of row names.\n");
  
  if (nrl>0)
   if (DEB & DEBJM)
    Rcpp::Rcout << "The passed matrix has row names for the " << nrl << " rows and they will be used.\n";
  
  if (mtype!=MTYPESYMMETRIC)
  {
   if (q4[1]!=R_NilValue)
   {
    colnames=q4[1];
    ncl=indextype(colnames.length());
   }
   else
    ncl=0;
     
   if ((ncl>0) && (ncl!=nc))
    Rcpp::stop("Strange Matrix object. The number of columns in the matrix differs from the length of the vector of column names.\n"); 

   if (ncl>0)
    if (DEB & DEBJM)
     Rcpp::Rcout << "The passed matrix has column names for the " << ncl << " columns and they will be used.\n";
  } 
 }
 
 switch (mtype)
 { 
  case MTYPEFULL:
  {
   FullMatrix<T> Mf(nr,nc);
   for (indextype r=0; r<nr; r++)
    for (indextype c=0; c<nc; c++)
     Mf.Set(r,c,M(r,c));
   
   if (comment!="")
    Mf.SetComment(comment);
   if (nrl>0)
    Mf.SetRowNames(rownames);
   if (ncl>0)
    Mf.SetColNames(colnames);
    
   Mf.WriteBin(fname); 
  } 
  break;
  case MTYPESPARSE:
  {
   SparseMatrix<T> Mf(nr,nc);
   for (indextype r=0; r<nr; r++)
    for (indextype c=0; c<nc; c++)
     Mf.Set(r,c,M(r,c));
   
   if (comment!="")
    Mf.SetComment(comment);
   if (nrl>0)
    Mf.SetRowNames(rownames);
   if (ncl>0)
    Mf.SetColNames(colnames);
    
   Mf.WriteBin(fname); 
  } 
  break;
  case MTYPESYMMETRIC:
  {
   SymmetricMatrix<T> Mf(nr);
   for (indextype r=0; r<nr; r++)
    for (indextype c=0; c<r+1; c++)
     Mf.Set(r,c,M(r,c));
   
   if (comment!="")
    Mf.SetComment(comment);
   if (nrl>0)
    Mf.SetRowNames(rownames);
   if (ncl>0)
    Mf.SetColNames(colnames); 
   
   Mf.WriteBin(fname); 
  }
  default: break;
 }
}

//' JWriteBin
//'
//' Writes a R matrix to a disk file as a binary matrix in the jmatrix format
//'
//'
//' Use this function cautiously. Differently to the functions to get one or more rows or columns from the binary file,
//' which book only the memory strictly needed for the vector/matrix and do not load all the binary file in memory,
//' this function books the full matrix in the requested data type and writes it later so with very big matrices
//' you might run out of memory.\cr
//' Type 'int' is really long int (8-bytes in most modern machines) so using 'int' or 'long' is equivalent.\cr
//' Type is coerced from double (the internal type of R matrices) to the requested type, which may provoke a loose of precision.\cr
//' If M is a named-R matrix, row and column names are written as metadata, too.\cr
//' Also, if you write as symmetric a matrix which is not such, only the lower-diagonal part will be written.
//' The rest of the data will be lost. In this case, if the matrix has row and column names, only row names are written.
//'
//' @param    M The R matrix to be written
//' @param    fname The name of the file to write
//' @param    dtype The data type of the matrix to be written: one of the strings 'short', 'int', 'long', 'float' or 'double'. Default: 'float'
//' @param    dmtype The matrix type: one of the strings 'full', 'sparse' or 'symmetric'. Default: 'full'
//' @param    comment A optional string with the comment to be added as metadata. Default: "" (empty string, no added comment)
//' @return   No return value, called for side effects (creates a file)
//' @examples
//' Rf <- matrix(runif(48),nrow=6)
//' rownames(Rf) <- c("A","B","C","D","E","F")
//' colnames(Rf) <- c("a","b","c","d","e","f","g","h")
//' tmpfile1=paste0(tempdir(),"/Rfullfloat.bin")
//' JWriteBin(Rf,tmpfile1,dtype="float",dmtype="full",comment="Full matrix of floats")
//' @export
// [[Rcpp::export]]
void JWriteBin(Rcpp::NumericMatrix M,std::string fname,std::string dtype="float",std::string dmtype="full",std::string comment="")
{
 unsigned char ctype=NOTYPE;
 if (dtype=="short")
  ctype=SSTYPE;
 if ((dtype=="int") || (dtype=="long"))
  ctype=SLTYPE;
 if (dtype=="float")
  ctype=FTYPE;
 if (dtype=="double")
  ctype=DTYPE;
 if (ctype==NOTYPE)
  Rcpp::stop("Allowed data types are only 'short', 'int', 'float' or 'double'.\n");
  
 unsigned char mtype=MTYPENOTYPE;
 if (dmtype=="full")
  mtype=MTYPEFULL;
 if (dmtype=="sparse")
  mtype=MTYPESPARSE;
 if (dmtype=="symmetric")
  mtype=MTYPESYMMETRIC;
 if (mtype==MTYPENOTYPE)
  Rcpp::stop("Parameter mtype must be one of the strings 'full', 'sparse' or 'symmetric'\n");
  
 switch (ctype)
 {
  case SSTYPE: { WriteRMatrixAsBin<short>(mtype,fname,M,comment); break;};
  case SLTYPE: { WriteRMatrixAsBin<long>(mtype,fname,M,comment); break;};
  case FTYPE:  { WriteRMatrixAsBin<float>(mtype,fname,M,comment); break;};
  case DTYPE:  { WriteRMatrixAsBin<double>(mtype,fname,M,comment); break;};
  default: break;
 }
}

