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

/************ Special functions to get the lower-diagonal part of a symmetric matrix as a vector *************/

// Auxiliary function to get the lower-diagonal part of a symmetric matrix

template <typename T>
void GSDiag(std::string fname,indextype nrows,Rcpp::NumericVector &v)
{
 T *data = new T [nrows];
  
 std::ifstream f(fname.c_str());
 
 // This is the beginning of row 1 in the binary symmetric data
 size_t offset=HEADER_SIZE+sizeof(T);
 f.seekg(offset,std::ios::beg);
 
 for (indextype r=1; r<nrows; r++)
 {
  // Here we read the r+1 values present in that row, including the (nr,nr) at the main diagonal (which will be normally 0 in a dissimilarity matrix)
  f.read((char *)data,(r+1)*sizeof(T));
  // but we copy all of them, except the last one, at the appropriate place of return array so that theay are ordered by column.
  for (indextype c=0; c<r; c++)
   v(c*(nrows-1)-((c*(c-1))/2)+r-c-1)=data[c];
 }
 
 f.close();
 delete[] data;
}

//' GetSubdiag
//' 
//' Takes a symmetric matrix and returns a vector with all its elements under the main diagonal (without those at the diagonal itself)
//' Done as an instrumental function to check the PAM in package cluster. To be removed in final version of the package.
//'
//' @param    fname The name of the file with the dissimilarity matrix in jmatrix binary format.
//' @return   The vector with the values under the main diagonal, sorted by columns (i.e.: m(2,1) .. m(n,1), m(3,2)..m(n,2),..., m(n-1,n))
//' @examples
//' Rns <- matrix(runif(49),nrow=7)
//' Rsym <- 0.5*(Rns+t(Rns))
//' rownames(Rsym) <- c("A","B","C","D","E","F","G")
//' colnames(Rsym) <- c("a","b","c","d","e","f","g")
//' tmpfile1=paste0(tempdir(),"/Rsymfloat.bin")
//' JWriteBin(Rsym,tmpfile1,dtype="float",dmtype="symmetric")
//' d<-GetSubdiag(tmpfile1)
//' Rsym
//' d
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector GetSubdiag(std::string fname)
{
 unsigned char mtype,ctype,endian,mdinfo;
 indextype nrows,ncols;
 
 MatrixType(fname,mtype,ctype,endian,mdinfo,nrows,ncols);
 if (mtype != MTYPESYMMETRIC)
  Rcpp::stop("This function admits only symmetric matrices.\n");
  
 Rcpp::NumericVector retv((nrows*(nrows-1))/2);
 switch (ctype)
 {
  case FTYPE:  { GSDiag<float>(fname,nrows,retv); break;};
  case DTYPE:  { GSDiag<double>(fname,nrows,retv); break;};
  case LDTYPE: { GSDiag<long double>(fname,nrows,retv); break;};
  default: Rcpp::stop("This function admits only matrices of float, double or long double.\n"); break;
 } 
 
 return(retv);
}
