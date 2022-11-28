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
#include <matmetadata.h>

extern unsigned char DEB;

/***********************************************************
 * 
 * Functions related with getting rows from the jmatrix binary format 
 *
 **********************************************************/
 
// Functions to read one or more columns. These are simpler, since jmatrix format is primarily row-oriented.

// Auxiliary functions to read one or many rows of the matrices stored in binary jmatrix format without reading the full matrix in memory

template <typename T>
void GetJustOneRowFromFull(std::string fname,indextype nr,indextype ncols,Rcpp::NumericVector &v)
{
 T *data = new T [ncols]; 
 
 std::ifstream f(fname.c_str());
 // Start of row nr is at the end of former rows, each of them having ncols elements
 f.seekg(HEADER_SIZE+nr*ncols*sizeof(T),std::ios::beg);
 // Here we simply read ncols elements
 f.read((char *)data,ncols*sizeof(T));
 f.close();
 
 for (size_t c=0; c<ncols; c++)
  v(c)=data[c];
  
 delete[] data;
}

// Auxiliary function to get from a full matrix the rows whose indexes are in vector nr. Rows are left consecutively in the passed matrix
// in the same order as in vector r, so if vector  is not ordered, resulting rows will be unordered, too.
template <typename T>
void GetManyRowsFromFull(std::string fname,std::vector<indextype> nr,indextype ncols,Rcpp::NumericMatrix &m)
{
 T *data = new T [ncols]; 
 
 std::ifstream f(fname.c_str());
 for (size_t t=0; t<nr.size(); t++)
 {
  // Start of row nr is at the end of former rows, each of them having ncols elements
  f.seekg(HEADER_SIZE+nr[t]*ncols*sizeof(T),std::ios::beg);
  // Here we simply read ncols elements
  f.read((char *)data,ncols*sizeof(T));
  // and put them in the matrix
  for (indextype c=0; c<ncols; c++)
   m(t,c)=data[c];
 }
 f.close();
 delete[] data;
}

template <typename T>
void GetJustOneRowFromSparse(std::string fname,indextype nr,indextype ncols,Rcpp::NumericVector &v)
{
 indextype ncr;
 
 std::ifstream f(fname.c_str());
 // Start of row nr is at the end of former rows, each of them having a different number of elements
 // we must go to the start of each row to find out how many...
 size_t offset=HEADER_SIZE;
 f.seekg(offset,std::ios::beg);
 f.read((char *)&ncr,sizeof(indextype));
 
 
 for (indextype r=0; r<nr; r++)
 {
  // Advance to jump over the ncr indextype indices, plus the own ncr value, plus the ncr values of data
  offset += ((ncr+1)*sizeof(indextype)+ncr*sizeof(T));
  f.seekg(offset,std::ios::beg);
  // At the beginning of row r: read how many element there are in it (ncr):
  f.read((char *)&ncr,sizeof(indextype));
 }

 // Clear the vector to be returned
 for (size_t c=0; c<ncols; c++)
  v(c)=0; 
  
 if (ncr!=0)
 {
  indextype *idata;
  T *data;
  // We are just at the beginning of the row we want to read, and we have already read its number of non-null entries (ncr)
  // Let's read its indices...
  idata = new indextype [ncr];
  f.read((char *)idata,ncr*sizeof(indextype));
  // ... and let's read the data
  data = new T [ncr];
  f.read((char *)data,ncr*sizeof(T));
 
  // Fill the appropriate places of the vector (those dictated by the indices in idata)
  for (size_t c=0; c<ncr; c++)
   v(idata[c])=data[c];
  
  delete[] data;
  delete[] idata;
 }
 
 f.close();
}

// WARNING: see if there is any problem. It seems to be corrected.
// Auxiliary function to get rows from a sparse matrix whose row indexes are in vector nr. Rows are left consecutively in the passed matrix
// in the same order as in vector nr, so if vector is not ordered, resulting rows will be unordered, too.
template <typename T>
void GetManyRowsFromSparse(std::string fname,std::vector<indextype> nr,indextype nrows,indextype ncols,Rcpp::NumericMatrix &m)
{
 std::vector<size_t> offsets(nrows);
 
 std::ifstream f(fname.c_str());

 indextype ncr;
 offsets[0]=HEADER_SIZE;
 for (size_t t=0;t<nrows;t++)
 {
  f.seekg(offsets[t],std::ios::beg);
  f.read((char *)&ncr,sizeof(indextype));
  // Advance to jump over the ncr indextype indices, plus the own ncr value, plus the ncr values of data
  if (t<nrows-1)
   offsets[t+1] = offsets[t]+((ncr+1)*sizeof(indextype)+ncr*sizeof(T));
 }
 
 // These are temporary arrays to store the indices and values of a row.
 // Since each row stores a different number of values, we book sufficient space for the biggest case, the number of columns
 // (this is the case of an absolutely unsparse row, no zeros in it)
 indextype *idata = new indextype [ncols];
 T *data = new T [ncols];
 
 for (size_t t=0; t<nr.size(); t++)
 {
  // Clear the corresponding row of the matrix to be returned
  for (size_t c=0; c<ncols; c++)
   m(t,c)=0;
  
  f.seekg(offsets[nr[t]],std::ios::beg);
  // At the beginning of row r: read how many element there are in it (ncr):
  f.read((char *)&ncr,sizeof(indextype));
  if (ncr!=0)
  {
   // Let's read its indices... No problem with size, ncr is always smaller than ncols
   f.read((char *)idata,ncr*sizeof(indextype));
   // ... and let's read the data
   f.read((char *)data,ncr*sizeof(T));
  }
  
  // Fill the appropriate places of the  (those dictated by the indices in idata)
  for (size_t c=0; c<ncr; c++)
   m(t,idata[c])=data[c];
 }
  
 delete[] data;
 delete[] idata;
 
 f.close();
}

template <typename T>
void GetJustOneRowFromSymmetric(std::string fname,indextype nr,indextype ncols,Rcpp::NumericVector &v)
{
 T *data = new T [ncols]; 
 
 std::ifstream f(fname.c_str());
 
 // This is the beginning of row nr in the binary symmetric data
 size_t offset=HEADER_SIZE+((nr*(nr+1))/2)*sizeof(T);
 f.seekg(offset,std::ios::beg);
 
 // Here we read the nr+1 values present in that row, including the (nr,nr) at the main diagonal (which will be normally 0 in a dissimilarity matrix)
 f.read((char *)data,(nr+1)*sizeof(T));
 
 // The rest of this row is not physically after; we must read the rest of the column that starts in the diagonal, down to the end.
 // It would be clearer to write r=nr+1; r<nrows; r++ but we haven't a variable nrows. But we don't need: symmetric matrices are square.
 offset=HEADER_SIZE+(nr+((nr+1)*(nr+2))/2)*sizeof(T);
 for (indextype r=nr+1; r<ncols; r++)
 {
  f.seekg(offset,std::ios::beg);
  // Here we read just one value...
  f.read((char *)&(data[r]),sizeof(T));
  // and advance to the next row, taking into account that each row has a different number of really stored columns, which is r+1
  offset += ((r+1)*sizeof(T));
 }
 f.close();
 
 for (size_t c=0; c<ncols; c++)
  v(c)=data[c];
  
 delete[] data;
}

// Auxiliary function to get from a symmetric matrix the rows whose indexes are in vector nr. Rows are left consecutively in the passed matrix
// in the same order as in vector r, so if vector  is not ordered, resulting rows will be unordered, too.
template <typename T>
void GetManyRowsFromSymmetric(std::string fname,std::vector<indextype> nr,indextype ncols,Rcpp::NumericMatrix &m)
{
 T *data = new T [ncols]; 
 
 std::ifstream f(fname.c_str());
 size_t offset;
 for (size_t t=0; t<nr.size(); t++)
 {
  // This is the beginning of row nr in the binary symmetric data
  offset=HEADER_SIZE+((nr[t]*(nr[t]+1))/2)*sizeof(T);
  f.seekg(offset,std::ios::beg);
  // Here we read the nr+1 values present in that row, including the (nr,nr) at the main diagonal (which will be normally 0 in a dissimilarity matrix)
  f.read((char *)data,(nr[t]+1)*sizeof(T));
 
  for (size_t c=0; c<nr[t]+1; c++)
   m(t,c)=data[c];
   
  // The rest of this row is not physically after; we must read the rest of the column that starts in the diagonal, down to the end.
  // It would be clearer to write r=nr+1; r<nrows; r++ but we haven't a variable nrows. But we don't need: symmetric matrices are square.
  offset=HEADER_SIZE+(nr[t]+((nr[t]+1)*(nr[t]+2))/2)*sizeof(T);
  for (indextype r=nr[t]+1; r<ncols; r++)
  {
   f.seekg(offset,std::ios::beg);
   // Here we read just one value...
   f.read((char *)&(data[r]),sizeof(T));
   // and advance to the next row, taking into account that each row has a different number of really stored columns, which is r+1
   offset += ((r+1)*sizeof(T));
  }
  
  for (size_t c=nr[t]+1; c<ncols; c++)
   m(t,c)=data[c];
  
 }
 
 f.close();
 delete[] data;
}

void OneRowFromAnything(std::string fname,unsigned char mtype,unsigned char ctype,indextype nr,indextype ncols,Rcpp::NumericVector &retv)
{
 if (mtype==MTYPEFULL)
 {
     switch (ctype)
     {
        case UCTYPE: { GetJustOneRowFromFull<unsigned char>(fname,nr,ncols,retv); break; };
        case SCTYPE: { GetJustOneRowFromFull<char>(fname,nr,ncols,retv); break;};
        case USTYPE: { GetJustOneRowFromFull<unsigned short>(fname,nr,ncols,retv); break;};
        case SSTYPE: { GetJustOneRowFromFull<short>(fname,nr,ncols,retv); break;};
        case UITYPE: { GetJustOneRowFromFull<unsigned int>(fname,nr,ncols,retv); break;};
        case SITYPE: { GetJustOneRowFromFull<int>(fname,nr,ncols,retv); break;};
        case ULTYPE: { GetJustOneRowFromFull<unsigned long>(fname,nr,ncols,retv); break;};
        case SLTYPE: { GetJustOneRowFromFull<long>(fname,nr,ncols,retv); break;};
        case FTYPE:  { GetJustOneRowFromFull<float>(fname,nr,ncols,retv); break;};
        case DTYPE:  { GetJustOneRowFromFull<double>(fname,nr,ncols,retv); break;};
        case LDTYPE: { GetJustOneRowFromFull<long double>(fname,nr,ncols,retv); break;};
        default: break;
    }
 }

 if (mtype==MTYPESPARSE)
 {
     switch (ctype)
     {
        case UCTYPE: { GetJustOneRowFromSparse<unsigned char>(fname,nr,ncols,retv); break; };
        case SCTYPE: { GetJustOneRowFromSparse<char>(fname,nr,ncols,retv); break; };
        case USTYPE: { GetJustOneRowFromSparse<unsigned short>(fname,nr,ncols,retv); break; };
        case SSTYPE: { GetJustOneRowFromSparse<short>(fname,nr,ncols,retv); break; };
        case UITYPE: { GetJustOneRowFromSparse<unsigned int>(fname,nr,ncols,retv); break; };
        case SITYPE: { GetJustOneRowFromSparse<int>(fname,nr,ncols,retv); break; };
        case ULTYPE: { GetJustOneRowFromSparse<unsigned long>(fname,nr,ncols,retv); break; };
        case SLTYPE: { GetJustOneRowFromSparse<long>(fname,nr,ncols,retv); break; };
        case FTYPE:  { GetJustOneRowFromSparse<float>(fname,nr,ncols,retv); break; }; 
        case DTYPE:  { GetJustOneRowFromSparse<double>(fname,nr,ncols,retv); break; };
        case LDTYPE: { GetJustOneRowFromSparse<long double>(fname,nr,ncols,retv); break; };
        default: break;
    }
 }
 
 if (mtype==MTYPESYMMETRIC)
 {
     switch (ctype)
     {
        case UCTYPE: { GetJustOneRowFromSymmetric<unsigned char>(fname,nr,ncols,retv); break; };
        case SCTYPE: { GetJustOneRowFromSymmetric<char>(fname,nr,ncols,retv); break;};
        case USTYPE: { GetJustOneRowFromSymmetric<unsigned short>(fname,nr,ncols,retv); break;};
        case SSTYPE: { GetJustOneRowFromSymmetric<short>(fname,nr,ncols,retv); break;};
        case UITYPE: { GetJustOneRowFromSymmetric<unsigned int>(fname,nr,ncols,retv); break;};
        case SITYPE: { GetJustOneRowFromSymmetric<int>(fname,nr,ncols,retv); break;};
        case ULTYPE: { GetJustOneRowFromSymmetric<unsigned long>(fname,nr,ncols,retv); break;};
        case SLTYPE: { GetJustOneRowFromSymmetric<long>(fname,nr,ncols,retv); break;};
        case FTYPE:  { GetJustOneRowFromSymmetric<float>(fname,nr,ncols,retv); break;}; 
        case DTYPE:  { GetJustOneRowFromSymmetric<double>(fname,nr,ncols,retv); break;};
        case LDTYPE: { GetJustOneRowFromSymmetric<long double>(fname,nr,ncols,retv); break;};
        default: break;
    }
 } 
}

void ManyRowsFromAnything(std::string fname,unsigned char mtype,unsigned char ctype,std::vector<indextype> nr,indextype nrows,indextype ncols,Rcpp::NumericMatrix &retm)
{
 if (mtype==MTYPEFULL)
 {
     switch (ctype)
     {
        case UCTYPE: { GetManyRowsFromFull<unsigned char>(fname,nr,ncols,retm); break; };
        case SCTYPE: { GetManyRowsFromFull<char>(fname,nr,ncols,retm); break;};
        case USTYPE: { GetManyRowsFromFull<unsigned short>(fname,nr,ncols,retm); break;};
        case SSTYPE: { GetManyRowsFromFull<short>(fname,nr,ncols,retm); break;};
        case UITYPE: { GetManyRowsFromFull<unsigned int>(fname,nr,ncols,retm); break;};
        case SITYPE: { GetManyRowsFromFull<int>(fname,nr,ncols,retm); break;};
        case ULTYPE: { GetManyRowsFromFull<unsigned long>(fname,nr,ncols,retm); break;};
        case SLTYPE: { GetManyRowsFromFull<long>(fname,nr,ncols,retm); break;};
        case FTYPE:  { GetManyRowsFromFull<float>(fname,nr,ncols,retm); break;};
        case DTYPE:  { GetManyRowsFromFull<double>(fname,nr,ncols,retm); break;};
        case LDTYPE: { GetManyRowsFromFull<long double>(fname,nr,ncols,retm); break;};
        default: break;
    }
 }

 if (mtype==MTYPESPARSE)
 {
     switch (ctype)
     {
        case UCTYPE: { GetManyRowsFromSparse<unsigned char>(fname,nr,nrows,ncols,retm); break; };
        case SCTYPE: { GetManyRowsFromSparse<char>(fname,nr,nrows,ncols,retm); break; };
        case USTYPE: { GetManyRowsFromSparse<unsigned short>(fname,nr,nrows,ncols,retm); break; };
        case SSTYPE: { GetManyRowsFromSparse<short>(fname,nr,nrows,ncols,retm); break; };
        case UITYPE: { GetManyRowsFromSparse<unsigned int>(fname,nr,nrows,ncols,retm); break; };
        case SITYPE: { GetManyRowsFromSparse<int>(fname,nr,nrows,ncols,retm); break; };
        case ULTYPE: { GetManyRowsFromSparse<unsigned long>(fname,nr,nrows,ncols,retm); break; };
        case SLTYPE: { GetManyRowsFromSparse<long>(fname,nr,nrows,ncols,retm); break; };
        case FTYPE:  { GetManyRowsFromSparse<float>(fname,nr,nrows,ncols,retm); break; }; 
        case DTYPE:  { GetManyRowsFromSparse<double>(fname,nr,nrows,ncols,retm); break; };
        case LDTYPE: { GetManyRowsFromSparse<long double>(fname,nr,nrows,ncols,retm); break; };
        default: break;
    }
 }
 
 if (mtype==MTYPESYMMETRIC)
 {
     switch (ctype)
     {
        case UCTYPE: { GetManyRowsFromSymmetric<unsigned char>(fname,nr,ncols,retm); break; };
        case SCTYPE: { GetManyRowsFromSymmetric<char>(fname,nr,ncols,retm); break;};
        case USTYPE: { GetManyRowsFromSymmetric<unsigned short>(fname,nr,ncols,retm); break;};
        case SSTYPE: { GetManyRowsFromSymmetric<short>(fname,nr,ncols,retm); break;};
        case UITYPE: { GetManyRowsFromSymmetric<unsigned int>(fname,nr,ncols,retm); break;};
        case SITYPE: { GetManyRowsFromSymmetric<int>(fname,nr,ncols,retm); break;};
        case ULTYPE: { GetManyRowsFromSymmetric<unsigned long>(fname,nr,ncols,retm); break;};
        case SLTYPE: { GetManyRowsFromSymmetric<long>(fname,nr,ncols,retm); break;};
        case FTYPE:  { GetManyRowsFromSymmetric<float>(fname,nr,ncols,retm); break;}; 
        case DTYPE:  { GetManyRowsFromSymmetric<double>(fname,nr,ncols,retm); break;};
        case LDTYPE: { GetManyRowsFromSymmetric<long double>(fname,nr,ncols,retm); break;};
        default: break;
    }
 } 
}

//' GetJRow
//'
//' Returns (as a R numeric vector) the requested row number from the matrix contained in a jmatrix binary file
//'
//' @param fname  String with the file name that contains the binary data.
//' @param nrow   The number of the row to be returned, in R-numbering (from 1)
//' @return       A numeric vector with the values of elements in the requested row
//' @examples
//' Rf <- matrix(runif(48),nrow=6)
//' rownames(Rf) <- c("A","B","C","D","E","F")
//' colnames(Rf) <- c("a","b","c","d","e","f","g","h")
//' tmpfile1=paste0(tempdir(),"/Rfullfloat.bin")
//' JWriteBin(Rf,tmpfile1,dtype="float",dmtype="full",comment="Full matrix of floats")
//' Rf[3,]
//' vf<-GetJRow(tmpfile1,3)
//' vf
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector GetJRow(std::string fname,int nrow)
{
 if (nrow<1)
  Rcpp::stop("Index in R-notation cannot be less than 1.\n");
  
 unsigned char mtype,ctype,endian,mdinfo;
 indextype nrows,ncols;
 indextype nr=indextype(nrow);
 
 MatrixType(fname,mtype,ctype,endian,mdinfo,nrows,ncols);
 if (nr>nrows)
  Rcpp::stop("Requested row is beyond the limit of the matrix.\n");
 
 nr--;  // This is to use C-index from now on.
 
 Rcpp::NumericVector retv(ncols);
 
 OneRowFromAnything(fname,mtype,ctype,nr,ncols,retv);
 
 if (mdinfo & COL_NAMES)
 {
  Rcpp::StringVector v=GetJColNames(fname);
  retv.names()=v;
 }
 
 return(retv);
}

//' GetJManyRows
//'
//' Returns (as a R numeric matrix) the rows with the requested row numbers from the matrix contained in a jmatrix binary file
//'
//' @param fname    String with the file name that contains the binary data.
//' @param extrows  A numeric vector with the indexes of the rows to be extracted, in R-numbering (from 1)
//' @return         A numeric matrix with the values of elements in the requested rows
//' @examples
//' Rf <- matrix(runif(48),nrow=6)
//' rownames(Rf) <- c("A","B","C","D","E","F")
//' colnames(Rf) <- c("a","b","c","d","e","f","g","h")
//' tmpfile1=paste0(tempdir(),"/Rfullfloat.bin")
//' JWriteBin(Rf,tmpfile1,dtype="float",dmtype="full",comment="Full matrix of floats")
//' Rf[c(1,4),]
//' vc<-GetJManyRows(tmpfile1,c(1,4))
//' vc
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix GetJManyRows(std::string fname,Rcpp::NumericVector extrows)
{  
 unsigned char mtype,ctype,endian,mdinfo;
 indextype nrows,ncols;
 
 MatrixType(fname,mtype,ctype,endian,mdinfo,nrows,ncols);
 std::vector<indextype> enr;
 for (int i=0; i<extrows.length(); i++)
 {
  if ( extrows(i)<1 || extrows(i)>nrows)
   Rcpp::stop("At least one of the requested rows is 0, or negative, or it is beyond the limit of the matrix.\n");
  enr.push_back(indextype(extrows(i)-1));     // This is to use C-index from now on.
 }
 
 Rcpp::NumericMatrix retm(enr.size(),ncols);
 
 ManyRowsFromAnything(fname,mtype,ctype,enr,nrows,ncols,retm);
 
 if (mdinfo & COL_NAMES)
 {
  Rcpp::StringVector v=GetJColNames(fname);
  colnames(retm)=v;
 }
 if (mdinfo & ROW_NAMES)
 {
  Rcpp::StringVector v=GetJRowNames(fname);
  Rcpp::StringVector vs(extrows.length());
  for (size_t t=0;t<size_t(extrows.length()); t++)
   vs(t)=v(extrows(t)-1);
  rownames(retm)=vs;
 }
 
 return(retm);
}

//' GetJRowByName
//'
//' Returns (as a R numeric vector) the requested named row from the matrix contained in a jmatrix binary file
//'
//' @param fname   String with the file name that contains the binary data.
//' @param rowname The name of the row to be returned. If the matrix has no row names, or the name is not found, an empty vector is returned
//' @return        A numeric vector with the values of elements in the requested row
//' @examples
//' Rf <- matrix(runif(48),nrow=6)
//' rownames(Rf) <- c("A","B","C","D","E","F")
//' colnames(Rf) <- c("a","b","c","d","e","f","g","h")
//' tmpfile1=paste0(tempdir(),"/Rfullfloat.bin")
//' JWriteBin(Rf,tmpfile1,dtype="float",dmtype="full",comment="Full matrix of floats")
//' Rf["C",]
//' vf<-GetJRowByName(tmpfile1,"C")
//' vf
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector GetJRowByName(std::string fname,std::string rowname)
{
 unsigned char mtype,ctype,endian,mdinfo;
 indextype nrows,ncols;
 
 MatrixType(fname,mtype,ctype,endian,mdinfo,nrows,ncols);
 if (!(mdinfo & ROW_NAMES))
 {
  Rcpp::warning("The matrix stored in that file has no row names as metadata. Returning empty vector.\n");
  return Rcpp::NumericVector();
 }
 
 Rcpp::StringVector v=GetJRowNames(fname);
 
 indextype nr=0;
 while (nr<(unsigned int)v.length() && v[nr]!=rowname)
  nr++;
 if (nr>=(unsigned int)v.length())
 {
  Rcpp::warning("Requested row name not found in the metadata. Returning empty vector.\n");
  return Rcpp::NumericVector();
 }
 
 Rcpp::NumericVector retv(ncols);
 
 OneRowFromAnything(fname,mtype,ctype,nr,ncols,retv);
 
 if (mdinfo & COL_NAMES)
 {
  Rcpp::StringVector cnames=GetJColNames(fname);
  retv.names()=cnames;
 }
 
 return(retv);
}

//' GetJManyRowsByNames
//'
//' Returns (as a R numeric matrix) the rows with the requested row names from the matrix contained in a jmatrix binary file
//'
//' @param fname        String with the file name that contains the binary data.
//' @param extrownames  A numeric vector with the names of the rows to be extracted. If the binary file has no row names, or _any_ of the row names is not present, an empty matrix is returned.
//' @return             A numeric matrix with the values of elements in the requested rows
//' @examples
//' Rf <- matrix(runif(48),nrow=6)
//' rownames(Rf) <- c("A","B","C","D","E","F")
//' colnames(Rf) <- c("a","b","c","d","e","f","g","h")
//' tmpfile1=paste0(tempdir(),"/Rfullfloat.bin")
//' JWriteBin(Rf,tmpfile1,dtype="float",dmtype="full",comment="Full matrix of floats")
//' Rf[c("A","C"),]
//' vf<-GetJManyRowsByNames(tmpfile1,c("A","C"))
//' vf
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix GetJManyRowsByNames(std::string fname,Rcpp::StringVector extrownames)
{  
 unsigned char mtype,ctype,endian,mdinfo;
 indextype nrows,ncols;
 
 MatrixType(fname,mtype,ctype,endian,mdinfo,nrows,ncols);
 if (!(mdinfo & ROW_NAMES))
 {
  Rcpp::warning("The matrix stored in that file has no row names as metadata. Returning empty matrix.\n");
  return Rcpp::NumericMatrix();
 }
 
 Rcpp::StringVector v=GetJRowNames(fname);
 
 std::vector<indextype> enr(extrownames.length());
 for (int i=0; i<extrownames.length(); i++)
 {
  indextype nr=0;
  while (nr<(unsigned int)v.length() && v[nr]!=extrownames(i))
   nr++;
  if (nr>=(unsigned int)v.length())
  {
   Rcpp::warning("At least one requested row name not found in the metadata. Returning empty matrix.\n");
   return Rcpp::NumericMatrix();
  }
  enr[i]=nr;
 }
 
 Rcpp::NumericMatrix retm(enr.size(),ncols);
 
 ManyRowsFromAnything(fname,mtype,ctype,enr,nrows,ncols,retm);
 
 if (mdinfo & COL_NAMES)
 {
  Rcpp::StringVector v=GetJColNames(fname);
  colnames(retm)=v;
 }

 rownames(retm)=extrownames;
 
 return(retm);
}

