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
 * Functions related with getting cols from the jmatrix binary format
 *
 **********************************************************/
 
// Functions to read one or more columns. These are more complex, since jmatrix format is primarily row-oriented.

// Auxiliary functions to read one or many columns of the matrices stored in binary jmatrix format without reading the full matrix in memory

template <typename T>
void GetJustOneColumnFromFull(std::string fname,indextype nc,indextype nrows,indextype ncols,Rcpp::NumericVector &v)
{
 T *data = new T [nrows]; 
 
 std::ifstream f(fname.c_str());
 // This places us at the nc column of first row (row 0)
 size_t offset=HEADER_SIZE+nc*sizeof(T);
 for (indextype r=0; r<nrows; r++)
 {
  f.seekg(offset,std::ios::beg);
  // We read one element
  f.read((char *)(&data[r]),sizeof(T));
  // This jumps exactly one row, up to just before the nc column of next row.
  // The number of bytes in one row is the number of columns multiplied by the size of one element.
  offset += (ncols*sizeof(T));
 }
 f.close();
 
 for (size_t r=0; r<nrows; r++)
  v(r)=data[r];
  
 delete[] data;
}

template <typename T>
void GetManyColumnsFromFull(std::string fname,std::vector<indextype> ncs,indextype nrows,indextype ncols,Rcpp::NumericMatrix &m)
{
 T data;
 
 std::ifstream f(fname.c_str());
 size_t offset;
 for (size_t t=0; t<ncs.size(); t++)
 {
  // This places us at the requested column of the first row
  offset=HEADER_SIZE+ncs[t]*sizeof(T);
  for (indextype r=0; r<nrows; r++)
  {
   f.seekg(offset,std::ios::beg);
   // We read one element
   f.read((char *)(&data),sizeof(T));
   m(r,t)=data;
   // This jumps exactly one row, up to just before the nc column of next row.
   // The number of bytes in one row is the number of columns multiplied by the size of one element.
   offset += (ncols*sizeof(T));
  }
  
 }
 f.close();
}

template <typename T>
void GetJustOneColumnFromSparse(std::string fname,indextype nc,indextype nrows,indextype ncols,Rcpp::NumericVector &v)
{
 T *data = new T [nrows];
 indextype *idata = new indextype [ncols];  // This is by excess. Most rows will not have ncols real columns, since the matrix is sparse.
 
 std::ifstream f(fname.c_str());
 
 // Start of row nr is at the end of former rows, each of them having a different number of elements
 // we must go to the start of each row to find out how many...
 size_t offset=HEADER_SIZE;
 indextype c,ncr;
 for (indextype r=0; r<nrows; r++)
 {
  f.seekg(offset,std::ios::beg);
  // At the beginning of row r: read how many element there are in it (ncr):
  f.read((char *)&ncr,sizeof(indextype));
  // Read the indices of this row
  f.read((char *)idata,ncr*sizeof(indextype));
  // See if the index of the column we are looking for is there (a while loop with premature exit is OK, indices are ordered)
  c=0;
  while ( (c<ncr) && (idata[c]<nc) )
   c++;
  if ((c>=ncr) || (idata[c]!=nc))
   data[r]=0;
  else
  {
   // offset is now the beginning of the current row. We jump the indices (and the number of non-null indices, which is the reason of '+1' in ncr+1)
   // and also jump the first c data, too.
   f.seekg(offset+((ncr+1)*sizeof(indextype))+(c*sizeof(T)),std::ios::beg);
   f.read((char *)(&data[r]),sizeof(T));
  }
  // We advance up to the beginning of next row. The size of this row consists on the nc indices, plus the number of non-null indices, plut the number of present values.
  offset += ((ncr+1)*sizeof(indextype)+(ncr*sizeof(T)));
 }
 
 f.close();
 
 for (size_t r=0; r<nrows; r++)
  v(r)=data[r];
  
 delete[] data;
 delete[] idata;
}

// Auxiliary function to get columns from a sparse matrix whose indexes are in vector nc. Columns are left consecutively in the passed matrix
// in the same order as in vector nc, so if vector is not ordered, resulting columns will be unordered, too.
template <typename T>
void GetManyColumnsFromSparse(std::string fname,std::vector<indextype> nc,indextype nrows,indextype ncols,Rcpp::NumericMatrix &m)
{
 // Differently to the function to get rows, here we need the offsets of absolutely all rows, since at least one element will be extracted from each of them
 std::vector<size_t> offsets(nrows,HEADER_SIZE);
 
 std::ifstream f(fname.c_str()); 
 size_t offset=HEADER_SIZE;
 
 indextype ncr;
 for (size_t t=0;t<nrows;t++)
 {
  offsets[t]=offset;
  f.seekg(offset,std::ios::beg);
  f.read((char *)&ncr,sizeof(indextype));
  offset += ((ncr+1)*sizeof(indextype)+ncr*sizeof(T));
 }
 
 // These are temporary arrays to store the indices and values of a row.
 // Since each row stores a different number of values, we book sufficient space for the biggest case, the number of columns
 // (this is the case of an absolutely unsparse row, no zeros in it)
 indextype *idata = new indextype [ncols];
 T *data = new T [ncols];
 
 // Here we run through all rows, getting from each of them the requested columns
 for (size_t t=0; t<nrows; t++)
 {
  f.seekg(offsets[t],std::ios::beg);
  // At the beginning of row r: read how many element there are in it (ncr):
  f.read((char *)&ncr,sizeof(indextype));
  // Let's read its indices... No problem with size, ncr is always smaller than ncols
  f.read((char *)idata,ncr*sizeof(indextype));
  // ... and let's read the data
  f.read((char *)data,ncr*sizeof(T));
   
  // Clear the corresponding column of the matrix to be returned
  for (size_t c=0; c<nc.size(); c++)
   m(t,c)=0;
  
  // Fill the appropriate places of the matrix (those dictated by the indices in idata)
  for (size_t c=0; c<nc.size(); c++)
  {
   size_t q=0;
   while ((q<ncr) && (idata[q]!=nc[c]))
    q++;
   if (q<ncr)
    m(t,c)=data[q];
  }
 }
 delete[] data;
 delete[] idata;
  
 f.close();
}


// GetJustOneColumnFromSymmetric is indeed the same as GetJustOneRowFromSymmetric, but it is easier to have it here by name
template <typename T>
void GetJustOneColumnFromSymmetric(std::string fname,indextype nr,indextype ncols,Rcpp::NumericVector &v)
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

// Again GetManyColumnsFromSymmetric is indeed the same as GetManyRowsFromSymmetric, but it is easier to have it here by name
// This explains the apparent inconsistency in the names of the variables.
// Auxiliary function to get from a symmetric matrix the columns whose indexes are in vector nr. Columns are left consecutively in the passed matrix
// in the same order as in vector r, so if vector is not ordered, resulting columns will be unordered, too.
template <typename T>
void GetManyColumnsFromSymmetric(std::string fname,std::vector<indextype> nr,indextype ncols,Rcpp::NumericMatrix &m)
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
   m(c,t)=data[c];                   // This is the difference w.r.t. the original, m is filled transposed
   
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
   m(c,t)=data[c];                  // This is the difference w.r.t. the original, m is filled transposed
 }
 
 f.close();
 
 delete[] data;
}

void OneColFromAnything(std::string fname,unsigned char mtype,unsigned char ctype,indextype nc,indextype nrows,indextype ncols,Rcpp::NumericVector &retv)
{ 
 if (mtype==MTYPEFULL)
 {
     switch (ctype)
     {
        case UCTYPE: { GetJustOneColumnFromFull<unsigned char>(fname,nc,nrows,ncols,retv); break; };
        case SCTYPE: { GetJustOneColumnFromFull<char>(fname,nc,nrows,ncols,retv); break;};
        case USTYPE: { GetJustOneColumnFromFull<unsigned short>(fname,nc,nrows,ncols,retv); break;};
        case SSTYPE: { GetJustOneColumnFromFull<short>(fname,nc,nrows,ncols,retv); break;};
        case UITYPE: { GetJustOneColumnFromFull<unsigned int>(fname,nc,nrows,ncols,retv); break;};
        case SITYPE: { GetJustOneColumnFromFull<int>(fname,nc,nrows,ncols,retv); break;};
        case ULTYPE: { GetJustOneColumnFromFull<unsigned long>(fname,nc,nrows,ncols,retv); break;};
        case SLTYPE: { GetJustOneColumnFromFull<long>(fname,nc,nrows,ncols,retv); break;};
        case FTYPE:  { GetJustOneColumnFromFull<float>(fname,nc,nrows,ncols,retv); break;};
        case DTYPE:  { GetJustOneColumnFromFull<double>(fname,nc,nrows,ncols,retv); break;};
        case LDTYPE: { GetJustOneColumnFromFull<long double>(fname,nc,nrows,ncols,retv); break;};
        default: break;
    }
 }

 if (mtype==MTYPESPARSE)
 {
     switch (ctype)
     {
        case UCTYPE: { GetJustOneColumnFromSparse<unsigned char>(fname,nc,nrows,ncols,retv); break; };
        case SCTYPE: { GetJustOneColumnFromSparse<char>(fname,nc,nrows,ncols,retv); break; };
        case USTYPE: { GetJustOneColumnFromSparse<unsigned short>(fname,nc,nrows,ncols,retv); break; };
        case SSTYPE: { GetJustOneColumnFromSparse<short>(fname,nc,nrows,ncols,retv); break; };
        case UITYPE: { GetJustOneColumnFromSparse<unsigned int>(fname,nc,nrows,ncols,retv); break; };
        case SITYPE: { GetJustOneColumnFromSparse<int>(fname,nc,nrows,ncols,retv); break; };
        case ULTYPE: { GetJustOneColumnFromSparse<unsigned long>(fname,nc,nrows,ncols,retv); break; };
        case SLTYPE: { GetJustOneColumnFromSparse<long>(fname,nc,nrows,ncols,retv); break; };
        case FTYPE:  { GetJustOneColumnFromSparse<float>(fname,nc,nrows,ncols,retv); break; }; 
        case DTYPE:  { GetJustOneColumnFromSparse<double>(fname,nc,nrows,ncols,retv); break; };
        case LDTYPE: { GetJustOneColumnFromSparse<long double>(fname,nc,nrows,ncols,retv); break; };
        default: break;
    }
 }
 
 if (mtype==MTYPESYMMETRIC)
 {
     switch (ctype)
     {
        case UCTYPE: { GetJustOneColumnFromSymmetric<unsigned char>(fname,nc,ncols,retv); break; };
        case SCTYPE: { GetJustOneColumnFromSymmetric<char>(fname,nc,ncols,retv); break;};
        case USTYPE: { GetJustOneColumnFromSymmetric<unsigned short>(fname,nc,ncols,retv); break;};
        case SSTYPE: { GetJustOneColumnFromSymmetric<short>(fname,nc,ncols,retv); break;};
        case UITYPE: { GetJustOneColumnFromSymmetric<unsigned int>(fname,nc,ncols,retv); break;};
        case SITYPE: { GetJustOneColumnFromSymmetric<int>(fname,nc,ncols,retv); break;};
        case ULTYPE: { GetJustOneColumnFromSymmetric<unsigned long>(fname,nc,ncols,retv); break;};
        case SLTYPE: { GetJustOneColumnFromSymmetric<long>(fname,nc,ncols,retv); break;};
        case FTYPE:  { GetJustOneColumnFromSymmetric<float>(fname,nc,ncols,retv); break;}; 
        case DTYPE:  { GetJustOneColumnFromSymmetric<double>(fname,nc,ncols,retv); break;};
        case LDTYPE: { GetJustOneColumnFromSymmetric<long double>(fname,nc,ncols,retv); break;};
        default: break;
    }
 } 
}

void ManyColumnsFromAnything(std::string fname,unsigned char mtype,unsigned char ctype,std::vector<indextype> nc,indextype nrows,indextype ncols,Rcpp::NumericMatrix &retm)
{
 if (mtype==MTYPEFULL)
 {
     switch (ctype)
     {
        case UCTYPE: { GetManyColumnsFromFull<unsigned char>(fname,nc,nrows,ncols,retm); break; };
        case SCTYPE: { GetManyColumnsFromFull<char>(fname,nc,nrows,ncols,retm); break;};
        case USTYPE: { GetManyColumnsFromFull<unsigned short>(fname,nc,nrows,ncols,retm); break;};
        case SSTYPE: { GetManyColumnsFromFull<short>(fname,nc,nrows,ncols,retm); break;};
        case UITYPE: { GetManyColumnsFromFull<unsigned int>(fname,nc,nrows,ncols,retm); break;};
        case SITYPE: { GetManyColumnsFromFull<int>(fname,nc,nrows,ncols,retm); break;};
        case ULTYPE: { GetManyColumnsFromFull<unsigned long>(fname,nc,nrows,ncols,retm); break;};
        case SLTYPE: { GetManyColumnsFromFull<long>(fname,nc,nrows,ncols,retm); break;};
        case FTYPE:  { GetManyColumnsFromFull<float>(fname,nc,nrows,ncols,retm); break;};
        case DTYPE:  { GetManyColumnsFromFull<double>(fname,nc,nrows,ncols,retm); break;};
        case LDTYPE: { GetManyColumnsFromFull<long double>(fname,nc,nrows,ncols,retm); break;};
        default: break;
    }
 }

 if (mtype==MTYPESPARSE)
 {
     switch (ctype)
     {
        case UCTYPE: { GetManyColumnsFromSparse<unsigned char>(fname,nc,nrows,ncols,retm); break; };
        case SCTYPE: { GetManyColumnsFromSparse<char>(fname,nc,nrows,ncols,retm); break; };
        case USTYPE: { GetManyColumnsFromSparse<unsigned short>(fname,nc,nrows,ncols,retm); break; };
        case SSTYPE: { GetManyColumnsFromSparse<short>(fname,nc,nrows,ncols,retm); break; };
        case UITYPE: { GetManyColumnsFromSparse<unsigned int>(fname,nc,nrows,ncols,retm); break; };
        case SITYPE: { GetManyColumnsFromSparse<int>(fname,nc,nrows,ncols,retm); break; };
        case ULTYPE: { GetManyColumnsFromSparse<unsigned long>(fname,nc,nrows,ncols,retm); break; };
        case SLTYPE: { GetManyColumnsFromSparse<long>(fname,nc,nrows,ncols,retm); break; };
        case FTYPE:  { GetManyColumnsFromSparse<float>(fname,nc,nrows,ncols,retm); break; }; 
        case DTYPE:  { GetManyColumnsFromSparse<double>(fname,nc,nrows,ncols,retm); break; };
        case LDTYPE: { GetManyColumnsFromSparse<long double>(fname,nc,nrows,ncols,retm); break; };
        default: break;
    }
 }
 
 if (mtype==MTYPESYMMETRIC)
 {
     switch (ctype)
     {
        case UCTYPE: { GetManyColumnsFromSymmetric<unsigned char>(fname,nc,nrows,retm); break; };
        case SCTYPE: { GetManyColumnsFromSymmetric<char>(fname,nc,nrows,retm); break;};
        case USTYPE: { GetManyColumnsFromSymmetric<unsigned short>(fname,nc,nrows,retm); break;};
        case SSTYPE: { GetManyColumnsFromSymmetric<short>(fname,nc,nrows,retm); break;};
        case UITYPE: { GetManyColumnsFromSymmetric<unsigned int>(fname,nc,nrows,retm); break;};
        case SITYPE: { GetManyColumnsFromSymmetric<int>(fname,nc,nrows,retm); break;};
        case ULTYPE: { GetManyColumnsFromSymmetric<unsigned long>(fname,nc,nrows,retm); break;};
        case SLTYPE: { GetManyColumnsFromSymmetric<long>(fname,nc,nrows,retm); break;};
        case FTYPE:  { GetManyColumnsFromSymmetric<float>(fname,nc,nrows,retm); break;}; 
        case DTYPE:  { GetManyColumnsFromSymmetric<double>(fname,nc,nrows,retm); break;};
        case LDTYPE: { GetManyColumnsFromSymmetric<long double>(fname,nc,nrows,retm); break;};
        default: break;
    }
 } 
}

//' GetJCol
//'
//' Returns (as a R numeric vector) the requested column number from the matrix contained in a jmatrix binary file
//'
//' @param fname  String with the file name that contains the binary data.
//' @param ncol   The number of the column to be returned, in R-numbering (from 1)
//' @return       A numeric vector with the values of elements in the requested column
//' @examples
//' Rf <- matrix(runif(48),nrow=6)
//' rownames(Rf) <- c("A","B","C","D","E","F")
//' colnames(Rf) <- c("a","b","c","d","e","f","g","h")
//' tmpfile1=paste0(tempdir(),"/Rfullfloat.bin")
//' JWriteBin(Rf,tmpfile1,dtype="float",dmtype="full",comment="Full matrix of floats")
//' Rf[,3]
//' vf<-GetJCol(tmpfile1,3)
//' vf
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector GetJCol(std::string fname,int ncol)
{
 if (ncol<1)
  Rcpp::stop("Index in R-notation cannot be less than 1.\n");
  
 unsigned char mtype,ctype,endian,mdinfo;
 indextype nrows,ncols;
 indextype nc=indextype(ncol);
 
 MatrixType(fname,mtype,ctype,endian,mdinfo,nrows,ncols);
 if (nc>ncols)
  Rcpp::stop("Requested column is beyond the limit of the matrix.\n");
 
 nc--;  // This is to use C-index from now on.
 
 Rcpp::NumericVector retv(nrows);
 
 OneColFromAnything(fname,mtype,ctype,nc,nrows,ncols,retv);

 if (mdinfo & ROW_NAMES)
 {
  Rcpp::StringVector vn=GetJRowNames(fname);
  retv.names()=vn;
 }
 
 return(retv);
}

//' GetJManyCols
//'
//' Returns (as a R numeric matrix) the columns with the requested column numbers from the matrix contained in a jmatrix binary file
//'
//' @param fname    String with the file name that contains the binary data.
//' @param extcols  A numeric vector with the indexes of the columns to be extracted, in R-numbering (from 1)
//' @return         A numeric matrix with the values of elements in the requested columns
//' @examples
//' Rf <- matrix(runif(48),nrow=6)
//' rownames(Rf) <- c("A","B","C","D","E","F")
//' colnames(Rf) <- c("a","b","c","d","e","f","g","h")
//' tmpfile1=paste0(tempdir(),"/Rfullfloat.bin")
//' JWriteBin(Rf,tmpfile1,dtype="float",dmtype="full",comment="Full matrix of floats")
//' vc<-GetJManyCols(tmpfile1,c(1,4))
//' vc
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix GetJManyCols(std::string fname,Rcpp::NumericVector extcols)
{
 unsigned char mtype,ctype,endian,mdinfo;
 indextype nrows,ncols;
 
 MatrixType(fname,mtype,ctype,endian,mdinfo,nrows,ncols);
 std::vector<indextype> enc;
 for (int i=0; i<extcols.length(); i++)
 {
  if ( extcols(i)<1 || extcols(i)>ncols)
   Rcpp::stop("At least one of the requested columns is 0, or negative, or it is beyond the limit of the matrix.\n");
  enc.push_back(indextype(extcols(i)-1));     // This is to use C-index from now on.
 }
 
 Rcpp::NumericMatrix retm(nrows,enc.size());
 
 ManyColumnsFromAnything(fname,mtype,ctype,enc,nrows,ncols,retm);
 
 if (mdinfo & ROW_NAMES)
 {
  Rcpp::StringVector v=GetJRowNames(fname);
  rownames(retm)=v;
 }
 if (mdinfo & COL_NAMES)
 {
  Rcpp::StringVector v=GetJColNames(fname);
  Rcpp::StringVector vs(extcols.length());
  for (size_t t=0;t<size_t(extcols.length()); t++)
   vs(t)=v(extcols(t)-1);
  colnames(retm)=vs;
 }
 
 return(retm);
}

//' GetJColByName
//'
//' Returns (as a R numeric vector) the requested named column from the matrix contained in a jmatrix binary file
//'
//' @param fname   String with the file name that contains the binary data.
//' @param colname The name of the column to be returned. If the matrix has no column names, or the name is not found, an empty vector is returned
//' @return        A numeric vector with the values of elements in the requested column
//' @examples
//' Rf <- matrix(runif(48),nrow=6)
//' rownames(Rf) <- c("A","B","C","D","E","F")
//' colnames(Rf) <- c("a","b","c","d","e","f","g","h")
//' tmpfile1=paste0(tempdir(),"/Rfullfloat.bin")
//' JWriteBin(Rf,tmpfile1,dtype="float",dmtype="full",comment="Full matrix of floats")
//' Rf[,"c"]
//' vf<-GetJColByName(tmpfile1,"c")
//' vf
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector GetJColByName(std::string fname,std::string colname)
{
 unsigned char mtype,ctype,endian,mdinfo;
 indextype nrows,ncols;
 
 MatrixType(fname,mtype,ctype,endian,mdinfo,nrows,ncols);
 if (!(mdinfo & COL_NAMES))
 {
  Rcpp::warning("The matrix stored in that file has no column names as metadata. Returning empty vector.\n");
  return Rcpp::NumericVector();
 }
 
 std::vector<std::string> rnames,cnames;
 InternalGetBinNames(fname,ROW_NAMES | COL_NAMES,rnames,cnames);
 
 indextype nc=0;
 while (nc<cnames.size() && cnames[nc]!=colname)
  nc++;
 if (nc>=cnames.size())
 {
  Rcpp::warning("Requested column name not found in the metadata. Returning empty vector.\n");
  return Rcpp::NumericVector();
 }
 
 Rcpp::NumericVector retv(nrows);

 OneColFromAnything(fname,mtype,ctype,nc,nrows,ncols,retv);
 
 if (mdinfo & ROW_NAMES)
  retv.names()=rnames;
 
 return(retv);
} 

//' GetJManyColsByNames
//'
//' Returns (as a R numeric matrix) the columns with the requested column names from the matrix contained in a jmatrix binary file
//'
//' @param fname        String with the file name that contains the binary data.
//' @param extcolnames  A numeric vector with the names of the columns to be extracted. If the binary file has no column names, or _any_ of the column names is not present, an empty matrix is returned.
//' @return             A numeric matrix with the values of elements in the requested columns
//' @examples
//' Rf <- matrix(runif(48),nrow=6)
//' rownames(Rf) <- c("A","B","C","D","E","F")
//' colnames(Rf) <- c("a","b","c","d","e","f","g","h")
//' tmpfile1=paste0(tempdir(),"/Rfullfloat.bin")
//' JWriteBin(Rf,tmpfile1,dtype="float",dmtype="full",comment="Full matrix of floats")
//' Rf[,c(1,4)]
//' vf<-GetJManyColsByNames(tmpfile1,c("a","d"))
//' vf
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix GetJManyColsByNames(std::string fname,Rcpp::StringVector extcolnames)
{
 unsigned char mtype,ctype,endian,mdinfo;
 indextype nrows,ncols;
 
 MatrixType(fname,mtype,ctype,endian,mdinfo,nrows,ncols);
 if (!(mdinfo & COL_NAMES))
 {
  Rcpp::warning("The matrix stored in that file has no column names as metadata. Returning empty matrix.\n");
  return Rcpp::NumericMatrix();
 }
 
 std::vector<std::string> rnames,cnames;
 InternalGetBinNames(fname,ROW_NAMES | COL_NAMES,rnames,cnames);
 
 std::vector<indextype> enc(extcolnames.length());
 for (int i=0; i<extcolnames.length(); i++)
 {
  indextype nc=0;
  while (nc<cnames.size() && cnames[nc]!=std::string(extcolnames(i)))
   nc++;
  if (nc>=cnames.size())
  {
   Rcpp::warning("At least one requested column name not found in the metadata. Returning empty matrix.\n");
   return Rcpp::NumericMatrix();
  }
  enc[i]=nc;
 }
 
 Rcpp::NumericMatrix retm(nrows,enc.size());
 
 ManyColumnsFromAnything(fname,mtype,ctype,enc,nrows,ncols,retm);
 
 if (mdinfo & ROW_NAMES)
 {
  Rcpp::StringVector v(rnames.size());
  for (indextype t=0;t<rnames.size();t++)
   v[t]=rnames[t];
  rownames(retm)=v;
 }
 
 colnames(retm)=extcolnames;
 
 return(retm);
}

