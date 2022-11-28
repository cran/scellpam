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

#include <Rcpp.h>
#include <sparsematrix.h>
#include <templatemacros.h>

extern unsigned char DEB;

template <typename T>
void sort_indexes_and_values(const std::vector<T> &v,std::vector<size_t> &idx,std::vector<indextype> &vdx)
{  
  for (size_t i = 0; i != idx.size(); ++i) idx[i] = i; 
  
  // sort indexes based on comparing values in v
  // using std::stable_sort instead of std::sort
  // to avoid unnecessary index re-orderings
  // when v contains elements of equal values 
  stable_sort(idx.begin(), idx.end(),
       [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

  for (size_t i = 0; i != idx.size(); ++i)
      vdx[i]=v[idx[i]];
}

template void sort_indexes_and_values(const std::vector<unsigned char> &v,std::vector<size_t> &idx,std::vector<indextype> &vdx);
template void sort_indexes_and_values(const std::vector<char> &v,std::vector<size_t> &idx,std::vector<indextype> &vdx);
template void sort_indexes_and_values(const std::vector<unsigned short> &v,std::vector<size_t> &idx,std::vector<indextype> &vdx);
template void sort_indexes_and_values(const std::vector<short> &v,std::vector<size_t> &idx,std::vector<indextype> &vdx);
template void sort_indexes_and_values(const std::vector<unsigned int> &v,std::vector<size_t> &idx,std::vector<indextype> &vdx);
template void sort_indexes_and_values(const std::vector<int> &v,std::vector<size_t> &idx,std::vector<indextype> &vdx);
template void sort_indexes_and_values(const std::vector<unsigned long> &v,std::vector<size_t> &idx,std::vector<indextype> &vdx);
template void sort_indexes_and_values(const std::vector<long> &v,std::vector<size_t> &idx,std::vector<indextype> &vdx);
template void sort_indexes_and_values(const std::vector<float> &v,std::vector<size_t> &idx,std::vector<indextype> &vdx);
template void sort_indexes_and_values(const std::vector<double> &v,std::vector<size_t> &idx,std::vector<indextype> &vdx);
template void sort_indexes_and_values(const std::vector<long double> &v,std::vector<size_t> &idx,std::vector<indextype> &vdx);

/****************************************
  TEMPLATED CONSTRUCTORS AND FUNCTIONS
*****************************************/

template <typename T>
SparseMatrix<T>::SparseMatrix() : JMatrix<T>(MTYPESPARSE)
{
 datacols.clear();
 data.clear();
}

TEMPLATES_CONST(SparseMatrix,)

////////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename T>
SparseMatrix<T>::SparseMatrix(indextype nrows,indextype ncols) : JMatrix<T>(MTYPESPARSE,nrows,ncols)
{
 std::vector<indextype> vc; vc.clear();
 std::vector<T> vt; vt.clear();
 
 for (indextype r=0;r<this->nr;r++)
 {
     datacols.push_back(vc);    
     data.push_back(vt);
 }
}

TEMPLATES_CONST(SparseMatrix,SINGLE_ARG(indextype nrows,indextype ncols))

////////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename T>
SparseMatrix<T>::SparseMatrix(const SparseMatrix<T>& other) : JMatrix<T>(other)
{
 if (this->nr != 0)
 {
  std::vector<indextype> vc; vc.clear();
  std::vector<T> vt; vt.clear();
 
  for (indextype r=0;r<this->nr;r++)
  {
     datacols.push_back(vc);    
     data.push_back(vt);
  }
  for (indextype r=0;r<this->nr;r++)
     for ( indextype c=0; c<other.datacols[r].size(); c++)
     {
       datacols[r].push_back(other.datacols[r][c]);  
       data[r].push_back(other.data[r][c]);
     }
 }
 else
 {
  datacols.clear();
  data.clear();   
 }
}

TEMPLATES_COPY_CONST(SparseMatrix)

////////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename T>
void SparseMatrix<T>::Resize(indextype newnr,indextype newnc)
{
 if (this->nr!=0)
  for (indextype r=0;r<this->nr;r++)
  {
    data[r].clear();
    datacols[r].clear();
  }
 data.clear();
 datacols.clear();
   
 ((JMatrix<T> *)this)->Resize(newnr,newnc);
 
 if (DEB & DEBJM)
       Rcpp::Rcout << "Sparse matrix resized to (" << this->nr << "," << this->nc << ")\n";
 
 std::vector<indextype> vc; vc.clear();
 std::vector<T> vt; vt.clear();
 
 for (indextype r=0;r<this->nr;r++)
 {
     datacols.push_back(vc);    
     data.push_back(vt);
 }
}

TEMPLATES_FUNC(void,SparseMatrix,Resize,SINGLE_ARG(indextype newnr,indextype newnc))

////////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename T>
SparseMatrix<T>::~SparseMatrix()
{
 if (this->nr!=0)
  for (indextype r=0;r<this->nr;r++)
  {
    data[r].clear();
    datacols[r].clear();
  }
 data.clear();
 datacols.clear();
}

TEMPLATES_DEFAULT_DEST(SparseMatrix)

////////////////////////////////////////////////////////////////////////////////////////////////////////

// Constructor to read from a binary file
template <typename T>
SparseMatrix<T>::SparseMatrix(std::string fname) : JMatrix<T>(fname,MTYPESPARSE)
{   
    std::vector<indextype> vc; vc.clear();
    std::vector<T> vt; vt.clear();
 
    for (indextype r=0;r<this->nr;r++)
    {
     datacols.push_back(vc);    
     data.push_back(vt);
    }
    
    indextype ncr;
    // These are booked by default to the maximum size of a row. Not all space will be used at each row.
    indextype *cvalues = new indextype [this->nc];
    T *values = new T [this->nc];
    
    // For each sparse row...
    for (indextype r=0;r<this->nr;r++)
    {
     // ...first, we read the number of non-zero entries
     this->ifile.read((char *)(&ncr),sizeof(indextype));
  
     // Then we read as many vales (column positions and real values) as needed, all at the same time, in arrays.
     this->ifile.read((char *)cvalues,ncr*sizeof(indextype));
     this->ifile.read((char *)values,ncr*sizeof(T));
     
     // and finally the arrays are stored as vectors
     for (indextype c=0;c<ncr;c++)
     {
         datacols[r].push_back(cvalues[c]);
         data[r].push_back(values[c]);
     }
    }
 
    delete[] cvalues;
    delete[] values;
    
    this->ReadMetadata();                  // This is exclusively used when reading from a binary file, not from a csv file
    
    this->ifile.close();
}

TEMPLATES_CONST(SparseMatrix,std::string fname)

////////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename T>
SparseMatrix<T>::SparseMatrix(std::string fname,TrMark) : JMatrix<T>(fname,MTYPESPARSE)
{
    // nr and nc have been read by the superclass contructor.
    // We need it, so we store them in orignr, orignc. But then we swap them
    indextype orignr,orignc;
    orignr=this->nr;
    orignc=this->nc;
    
    this->nr=orignc;
    this->nc=orignr;
    
    std::vector<indextype> vc; vc.clear();
    std::vector<T> vt; vt.clear();
 
    for (indextype r=0;r<this->nr;r++)
    {
     datacols.push_back(vc);    
     data.push_back(vt);
    }
    
    indextype ncr;
    // These are booked by default to the maximum size of a row. Not all space will be used at each row.
    indextype *cvalues = new indextype [orignc];
    T *values = new T [orignc];
    
    // For each sparse row...
    for (indextype r=0;r<orignr;r++)
    {
     // ...first, we read the number of non-zero entries
     this->ifile.read((char *)(&ncr),sizeof(indextype));
  
     // Then we read as many vales (column positions and real values) as needed, all at the same time, in arrays.
     this->ifile.read((char *)cvalues,ncr*sizeof(indextype));
     this->ifile.read((char *)values,ncr*sizeof(T));
     
     // and the arrays are stored as vectors in this new matrix, but transposed:
     // all are in row r, and in columns cvalues[0], cvalues[1], etc. so they will go to row cvalues[..], column r.
     for (indextype c=0;c<ncr;c++)
     {
         datacols[cvalues[c]].push_back(r);
         data[cvalues[c]].push_back(values[c]);
     }
    }
    
    delete[] cvalues;
    delete[] values;
    
    this->ReadMetadata();    // This is exclusively used when reading from a binary file, not from a csv file        
    
    this->ifile.close();
    
    // But this is not enough: each datum at (r,c) has been stored in row c, but not ordered, as expected. So we have to order each new row.
    // This time r runs up to the number of rows of the matrix that is being created, nr, not up to orignr as before
    for (indextype r=0;r<this->nr;r++)
    {
        std::vector<size_t>     idx(datacols[r].size());
        std::vector<indextype>  idv(datacols[r].size());
        // idv are the column-values, correctly ordered.
        // idx are their indexes, which are needed to reorder the array of values...
        sort_indexes_and_values(datacols[r],idx,idv);
        datacols[r].clear();
        datacols[r]=idv;
        std::vector<T> values;
        for (indextype c=0;c<idx.size();c++)
            values.push_back(data[r][idx[c]]);
        data[r].clear();
        data[r]=values;
    }
}

TEMPLATES_CONST(SparseMatrix,SINGLE_ARG(std::string fname,TrMark))

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename T>
SparseMatrix<T>& SparseMatrix<T>::operator=(const SparseMatrix<T>& other)
{
 if (this->nr!=0)
 {
  for (indextype r=0;r<this->nr;r++)
  {
   data[r].clear(); 
   datacols[r].clear();
  }
  data.clear();
  datacols.clear();
 }
 
 ((JMatrix<T> *)this)->operator=((const JMatrix<T> &)other);
 
 std::vector<indextype> vc; vc.clear();
 std::vector<T> vt; vt.clear();
 
 for (indextype r=0;r<this->nr;r++)
 {
     datacols.push_back(vc);    
     data.push_back(vt);
 }
 
 
 for (indextype r=0;r<this->nr;r++)
     for ( indextype c=0; c<other.datacols[r].size(); c++)
     {
       datacols[r].push_back(other.datacols[r][c]);  
       data[r].push_back(other.data[r][c]);
     }  
 
 return *this;
}

TEMPLATES_OPERATOR(SparseMatrix,=)

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename T>
SparseMatrix<T>& SparseMatrix<T>::operator!=(const SparseMatrix<T>& other)
{
 if (this->nr!=0)
 {
  if (DEB & DEBJM)
      Rcpp::Rcout << "Cleaning old matrix before assignment...\n";
  for (indextype r=0;r<this->nr;r++)
  {
   data[r].clear(); 
   datacols[r].clear();
  }
  data.clear();
  datacols.clear();
 }
 
 ((JMatrix<T> *)this)->operator!=((const JMatrix<T> &)other);
 // Here number of rows and columns has been swapped
 
 if (DEB & DEBJM)
 {
     indextype oldr,oldc;
     oldr=((JMatrix<T> *)&other)->GetNRows();
     oldc=((JMatrix<T> *)&other)->GetNCols();
     Rcpp::Rcout << "Transposing matrix of (" << oldr << "x" << oldc << ") to a matrix of (" << this->nr << "x" << this->nc << ")\n";
 }
 std::vector<indextype> vc; vc.clear();
 std::vector<T> vt; vt.clear();
 
 for (indextype r=0;r<this->nr;r++)
 {
     datacols.push_back(vc);    
     data.push_back(vt);
 }
  
 T v;
 for (indextype r=0;r<this->nr;r++)
 {
  for ( indextype c=0; c<this->nc; c++)
  {
    v=other.Get(c,r);
    if (v!=T(0))
    {
     datacols[r].push_back(c);
     data[r].push_back(v);   
    }
  }
 }
 return *this;
}

TEMPLATES_OPERATOR(SparseMatrix,!=)

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Constructor to read from a csv file
template <typename T>
SparseMatrix<T>::SparseMatrix(std::string fname,unsigned char vtype,char csep) : JMatrix<T>(fname,MTYPESPARSE,vtype,csep)
{    
    std::string line;
    // This is just to know number of rows
    this->nr=0;
    while (!this->ifile.eof())
    {
        this->ifile >> line;
        if (!this->ifile.eof())
            this->nr++;
    }
    if (DEB & DEBJM)
    {
        Rcpp::Rcout << this->nr << " lines (excluding header) in file " << fname << std::endl;
        Rcpp::Rcout << "Data will be read from each line and stored as ";
        switch (vtype)
        {
         case ULTYPE: Rcpp::Rcout << "unsigned 32-bit integers.\n"; break;
         case FTYPE: Rcpp::Rcout << "float values.\n"; break;
         case DTYPE: Rcpp::Rcout << "double values.\n"; break;
         default: Rcpp::Rcout << "unknown type values??? (Is this an error?).\n"; break;
        }
    }
    // Reposition pointer at the begin and re-read first (header) line
    // and no, seekg does not work (possibly, because we had reached the eof ???)
    this->ifile.close();
    this->ifile.open(fname.c_str());
    size_t p=0;
    this->ifile >> line;
    // No need to process first line here, it was done at the parent's class constructor
    
    T *data_with_zeros = new T [this->nc];
    
    std::vector<indextype> datacolsofrow;
    std::vector<T>dataofrow;
    
    while (!this->ifile.eof())
    {
        this->ifile >> line;
        if (!this->ifile.eof())
        {
          if (!this->ProcessDataLineCsv(line,csep,data_with_zeros))
          {
              std::ostringstream errst;
              errst << "Format error reading line " << p << " of file " << fname << ".\n";
              Rcpp::stop(errst.str());
          }
          datacolsofrow.clear();
          dataofrow.clear();
          for (indextype t=0; t<this->nc; t++)
              if (data_with_zeros[t]!=0)
              {
                  datacolsofrow.push_back(t);
                  dataofrow.push_back(data_with_zeros[t]);
              }
          datacols.push_back(datacolsofrow);
          data.push_back(dataofrow);
          
          p++;
        }
    }
    delete [] data_with_zeros;
     
    if (DEB & DEBJM)
    {
        Rcpp::Rcout << "Read " << p << " data lines of file " << fname;
        if (p==this->nr)
            Rcpp::Rcout << ", as expected.\n";
        else
            Rcpp::Rcout << " instead of " << this->nr << ".\n";
    }
    
    // No call to ReadMetadata must be done here, since these data are NOT binary. The column names were read by the parent's class constructor
    // and the row names are stored as they are read by the former call to ProcessDataLine
    
    this->ifile.close();  
}

TEMPLATES_CONST(SparseMatrix,SINGLE_ARG(std::string fname,unsigned char vtype,char csep))

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename T>
T SparseMatrix<T>::Get(indextype r,indextype c) const
{ 
#ifdef WITH_CHECKS_MATRIXSP
    if ((r>=this->nr) || (c>=this->nc))
    {
        std::ostringstream errst;
        errst << "Runtime error in SparseMatrix<T>::Get: at least one index (" << r << " or " << c << ") out of bounds.\n";
        errst << "This matrix was of dimension (" << this->nr << " x " << this->nc << ")\n";
        Rcpp::stop(errst.str());
    }
#endif
    // If this row has nothing in it or the first column in which there is something is beyond our request, the content is zero.
    if ((datacols[r].size()==0) || (datacols[r][0]>c))  // This assumes short-cut evaluation of OR...
        return T(0);

    // Otherwise, let's look for our column, if it exists.
    size_t m;
    size_t mid;
    size_t low=0;
    size_t high=datacols[r].size()-1;
    while (low <= high)
    {
     mid = low + (high-low)/2;
     
     if (datacols[r][mid] == c)
     {
         m = mid;
         break;
     }
     if (datacols[r][mid] < c)
         low=mid+1;
     else
         high=mid-1;
    }
    if (low>high)
        m = mid;
    
    return (datacols[r][m]==c) ? data[r][m] : T(0);
   
}

TEMPLATES_FUNCRCONST(SparseMatrix,Get,SINGLE_ARG(indextype r,indextype c))

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename T>
void SparseMatrix<T>::Set(indextype r,indextype c,T v)
{
#ifdef WITH_CHECKS_MATRIXSP
    if ((r>=this->nr) || (c>=this->nc))
    {
    	std::ostringstream errst;
        errst << "Runtime error in SparseMatrix<T>::Set: at least one index (" << r << " or " << c << ") out of bounds.\n";
        errst << "This matrix was of dimension (" << this->nr << " x " << this->nc << ")\n";
        Rcpp::stop(errst.str());
    }
#endif
    // First of all: zeros are not added.
    if (v==T(0))
     return;
     
    // Particular case 1: If this row has nothing in it, just add column and value
    if (datacols[r].size()==0)
    {
       datacols[r].push_back(c);
       data[r].push_back(v);
       return;
    } 
    
    // Particular case 2: It the first column in which there is something is beyond our request, the content must be set at the beginning
    if (datacols[r][0]>c)
    {
     // Yes, insert at head is begin()+1, not begin...
     datacols[r].insert(datacols[r].begin()+1,c);
     data[r].insert(data[r].begin()+1,v);
    }
    else
    {
     size_t mid;
     size_t low=0;
     size_t high=datacols[r].size()-1;
     while (low <= high)
     {
      mid = low + (high-low)/2;
     
      if (datacols[r][mid] == c)
      {
         data[r][mid]=v;
         return;
      }
      if (datacols[r][mid] < c)
          low=mid+1;
      else
          high=mid-1;
     }
     
     // If we are here, the column did not exist. m will contain the place to be inserted. Rememeber that insert after mid is at begin()+mid+1
     datacols[r].insert(datacols[r].begin()+mid+1,c);
     data[r].insert(data[r].begin()+mid+1,v);
    } 
    return;
}

TEMPLATES_SETFUNC(void,SparseMatrix,Set,SINGLE_ARG(indextype r,indextype c),v)

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename T>
void SparseMatrix<T>::SetRow(indextype r,std::vector<indextype> vc,std::vector<T> v)
{
#ifdef WITH_CHECKS_MATRIXSP
    if ((r>=this->nr) || (vc.size()>=this->nc))
    {
        std::ostringstream errst;
        errst << "Runtime error in SparseMatrix<T>::SetRow: either the row index " << r << " or the size of vc, " << vc.size() << " is/are out of bounds.\n";
        errst << "This matrix was of dimension (" << this->nr << " x " << this->nc << ")\n";
        Rcpp::stop(errst.str());
    }
#endif
    datacols[r].clear();
    datacols[r]=vc;
    data[r].clear();
    data[r]=v;
}

TEMPLATES_SETFUNCVEC(void,SparseMatrix,SetRow,SINGLE_ARG(indextype r,std::vector<indextype> vc),v)

//////////////////////////////////////////////////////////////////////////////////////////
template <typename T>
void SparseMatrix<T>::GetRow(indextype r,T *v)
{
#ifdef WITH_CHECKS_MATRIXSP
    if (r>=this->nr)
    {
    	std::ostringstream errst;
        errst << "Runtime error in SparseMatrix<T>::GetRow: the row index " << r << " is out of bounds.\n";
        errst << "This matrix was of dimension (" << this->nr << " x " << this->nc << ")\n";
        Rcpp::stop(errst.str());
    }
#endif
 // Fill the positions in v which are not zero.
 for (indextype c=0;c<data[r].size();c++)
     v[datacols[r][c]]=data[r][c];
}

TEMPLATES_SETFUNC(void,SparseMatrix,GetRow,indextype r,*v)

///////////////////////////////////////////////////////////////////////////////////////////////////

template <typename T>
void SparseMatrix<T>::GetSparseRow(indextype r,unsigned char *m,unsigned char s,T *v)
{
#ifdef WITH_CHECKS_MATRIXSP
    if (r>=this->nr)
    {
        std::ostringstream errst;
        errst << "Runtime error in SparseMatrix<T>::GetSparseRow: the row index " << r << " is out of bounds.\n";
        errst << "This matrix was of dimension (" << this->nr << " x " << this->nc << ")\n";
        Rcpp::stop(errst.str());
    }
#endif
  // Fill the positions in v which are not zero and also sum the value s to those positions in array m
  for (indextype c=0;c<data[r].size();c++)
  {
     v[datacols[r][c]]=data[r][c];  
     m[datacols[r][c]] |= s;
  }
}

TEMPLATES_SETFUNC(void,SparseMatrix,GetSparseRow,SINGLE_ARG(indextype r,unsigned char *m,unsigned char s),*v)

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename T>
void SparseMatrix<T>::GetMarksOfSparseRow(indextype r,unsigned char *m,unsigned char s)
{
#ifdef WITH_CHECKS_MATRIXSP
    if (r>=this->nr)
    {
    	std::ostringstream errst;
        errst << "Runtime error in SparseMatrix<T>::GetSparseRow: the row index " << r << " is out of bounds.\n";
        errst << "This matrix was of dimension (" << this->nr << " x " << this->nc << ")\n";
        Rcpp::stop(errst.str());
    }
#endif
  for (indextype c=0;c<data[r].size();c++)
     m[datacols[r][c]] |= s;
}

TEMPLATES_FUNC(void,SparseMatrix,GetMarksOfSparseRow,SINGLE_ARG(indextype r,unsigned char *m,unsigned char s))

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename T>
void SparseMatrix<T>::SelfRowNorm(std::string ctype)
{ 
 if (DEB & DEBJM)
  Rcpp::Rcout << "Normalizing... ";

 if ((ctype=="log1") || (ctype=="log1n"))
 {
  for (indextype r=0;r<this->nr;r++)
   for (indextype k=0;k<datacols[r].size();k++)
    data[r][k] = log2(data[r][k]+1.0);
 }
 
 if (ctype=="log1")
 {
  if (DEB & DEBJM)
   Rcpp::Rcout << "done!\n";
  return;
 }
 
 T sum;
 for (indextype r=0;r<this->nr;r++)
 {
  sum=T(0);
  for (indextype k=0;k<datacols[r].size();k++)
   sum+=data[r][k];
 
 if (sum!=T(0))
  for (indextype k=0;k<datacols[r].size();k++)
    data[r][k] /= sum;
 }
 if (DEB & DEBJM)
   Rcpp::Rcout << "done!\n";
}

TEMPLATES_FUNC(void,SparseMatrix,SelfRowNorm,std::string ctype)

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename T>
void SparseMatrix<T>::SelfColNorm(std::string ctype)
{ 
 if (DEB & DEBJM)
  Rcpp::Rcout << "Normalizing... ";

 if ((ctype=="log1") || (ctype=="log1n"))
 {
  for (indextype r=0;r<this->nr;r++)
   for (indextype k=0;k<datacols[r].size();k++)
    data[r][k] = log2(data[r][k]+1.0);
 }
 
 if (ctype=="log1")
 {
  if (DEB & DEBJM)
   Rcpp::Rcout << "done!\n";
  return;
 }
 
 T *sums = new T [this->nc];
 for (indextype c=0; c<this->nc; c++)
  sums[c]=T(0);
  
 for (indextype r=0; r<this->nr; r++)
  for (indextype k=0; k<datacols[r].size(); k++)
   sums[datacols[r][k]] += data[r][k];
 
 for (indextype r=0; r<this->nr; r++)
  for (indextype k=0; k<datacols[r].size(); k++)
   if (datacols[r][k]!=T(0))
    data[r][k] /= sums[datacols[r][k]];
 
 delete[] sums;
   
 if (DEB & DEBJM)
   Rcpp::Rcout << "done!\n";
}

TEMPLATES_FUNC(void,SparseMatrix,SelfColNorm,std::string ctype)

////////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename T>
void SparseMatrix<T>::WriteBin(std::string fname)
{
    ((JMatrix<T> *)this)->WriteBin(fname,MTYPESPARSE);
    
    if (DEB & DEBJM)
    {
     Rcpp::Rcout << "Writing binary matrix " << fname << " of (" << this->nr << "x" << this->nc << ")\n";
     Rcpp::Rcout.flush();
    }
    
    indextype ncr;
    for (indextype r=0;r<this->nr;r++)
    {
        ncr=datacols[r].size();
        this->ofile.write((const char *)(&ncr),sizeof(indextype));
        for (unsigned c=0;c<ncr;c++)
            this->ofile.write((const char *)(&datacols[r][c]),sizeof(indextype));
        for (unsigned c=0;c<ncr;c++)
            this->ofile.write((const char *)(&data[r][c]),sizeof(T));
    }
    
    unsigned long long endofbindata = this->ofile.tellp();
    
    if (DEB & DEBJM)
     Rcpp::Rcout << "End of block of binary data at offset " << endofbindata << "\n";

    this->WriteMetadata();              // Here we must write the metadata at the end of the binary contents of the matrix
    
    this->ofile.write((const char *)&endofbindata,sizeof(unsigned long long));  // This writes the point where binary data ends at the end of the file
    
    this->ofile.close();
}

TEMPLATES_FUNC(void,SparseMatrix,WriteBin,std::string fname)

/////////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename T>
void SparseMatrix<T>::WriteCsv(std::string fname,char csep,bool withquotes)
{
    ((JMatrix<T> *)this)->WriteCsv(fname,csep,withquotes);
    
    bool with_headers=false;
    size_t nch=this->colnames.size();
    size_t nrh=this->rownames.size();
    
    if (nch>0 && nrh>0)
    {
     if (nch!=this->nc || nrh!=this->nr)
      Rcpp::warning("Different size of headers and matrix, either in rows or in columns. Headers will not be written in the .csv file.\n");
     with_headers=true;
    }
   
    for (indextype r=0;r<this->nr;r++)
    {
        if (with_headers)
            this->ofile << FixQuotes(this->rownames[r],withquotes) << csep;
           
        for (indextype c=0;c<this->nc-1;c++)
            this->ofile << Get(r,c) << csep;
        this->ofile << Get(r,this->nc-1) << std::endl;
    }
    this->ofile.close();
}

TEMPLATES_FUNC(void,SparseMatrix,WriteCsv,SINGLE_ARG(std::string fname,char csep,bool withquotes))

/////////////////////////////////////////////////////////////////////////////////////////

template <typename T>
float SparseMatrix<T>::GetUsedMemoryMB()
{
    unsigned long long num_elem=0;
    for (indextype r=0;r<this->nr;r++)
        num_elem += (unsigned long long)datacols[r].size();
        
    Rcpp::Rcout << num_elem << " elements, half of " << sizeof(T) << " bytes and half of " << sizeof(indextype) << " bytes each, with accounts for ";
    float ret = float(num_elem)*float(sizeof(indextype)+sizeof(T));
    ret += float(datacols.size()*sizeof(indextype));
    ret /= (1024.0*1024.0);
    return ret;
}

TEMPLATES_FUNC(float,SparseMatrix,GetUsedMemoryMB,)

