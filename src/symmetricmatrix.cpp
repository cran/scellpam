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

#include <symmetricmatrix.h>


extern unsigned char DEB;

/****************************************
  TEMPLATED CONSTRUCTORS AND FUNCTIONS
*****************************************/

template <typename T>
SymmetricMatrix<T>::SymmetricMatrix() : JMatrix<T>(MTYPESYMMETRIC)
{
 data.clear();
}

TEMPLATES_CONST(SymmetricMatrix,)

//////////////////////////////////////////////////////////////////////////////////////////////

template <typename T>
SymmetricMatrix<T>::SymmetricMatrix(indextype nrows) : JMatrix<T>(MTYPESYMMETRIC,nrows,nrows)
{
 data.resize(this->nr);
 for (indextype r=0;r<this->nr;r++)
 {
     data[r].resize(r+1);
     data[r].assign(r+1,T(0));
 
 }
}

TEMPLATES_CONST(SymmetricMatrix,indextype rows)

//////////////////////////////////////////////////////////////////////////////////////////////

template <typename T>
SymmetricMatrix<T>::SymmetricMatrix(const SymmetricMatrix& other) : JMatrix<T>(other)
{
 data.resize(this->nr);
 for (indextype r=0;r<this->nr;r++)
 {
     data[r].resize(r+1);
     copy(other.data[r].begin(),other.data[r].end(),data[r].begin());
  
 }
}

TEMPLATES_COPY_CONST(SymmetricMatrix)

//////////////////////////////////////////////////////////////////////////////////////////////

template <typename T>
void SymmetricMatrix<T>::Resize(indextype newnr)
{
   if (data.size()!=0)
   {
    for (indextype r=0;r<data.size();r++)
     data[r].clear();
   }
   
   ((JMatrix<T> *)this)->Resize(newnr,newnr);
   
   if (DEB & DEBJM)
       Rcpp::Rcout << "Symmetric matrix resized to (" << this->nr << "," << this->nc << ")\n";
   
   data.resize(this->nr);
   for (indextype r=0;r<this->nr;r++)
   {
     data[r].resize(r+1);
     for (indextype c=0;c<=r;c++)
         data[r][c]=T(0);
   }
}

TEMPLATES_FUNC(void,SymmetricMatrix,Resize,indextype newnr)

//////////////////////////////////////////////////////////////////////////////////////////////

template <typename T>
SymmetricMatrix<T>::~SymmetricMatrix()
{
 if (data.size()!=0)
 {
  for (indextype r=0;r<data.size();r++)
   data[r].clear();
 }
}

TEMPLATES_DEFAULT_DEST(SymmetricMatrix)

//////////////////////////////////////////////////////////////////////////////////////////////

// Constructor to read from a binary file
template <typename T>
SymmetricMatrix<T>::SymmetricMatrix(std::string fname) : JMatrix<T>(fname,MTYPESYMMETRIC)
{
    try
    {
     data.resize(this->nr);
    }
    catch (int e)
    {
        Rcpp::stop("Exception allocating data.\n");
    }
    for (indextype r=0;r<this->nr;r++)
    {
        try
        {
         data[r].resize(r+1);
        }
        catch (int e)
        {
         std::ostringstream errst;
         errst << "Exception allocating data[" << r << "].\n";
         Rcpp::stop(errst.str());
        }
    }
    
    T *ddata = new T [this->nr];
    
    for (indextype r=0;r<this->nr;r++)
    {
        this->ifile.read((char *)ddata,(r+1)*sizeof(T));
        for (indextype c=0;c<=r;c++)
            data[r][c]=ddata[c];
    }
    
    delete[] ddata;
    
    this->ReadMetadata();                  // This is exclusively used when reading from a binary file, not from a csv file
      
    this->ifile.close();
    
    if (DEB & DEBJM)
     Rcpp::Rcout << "Read symmetric matrix with size (" << this->nr << "," << this->nc << ")\n";
     
}

TEMPLATES_CONST(SymmetricMatrix,std::string fname)

//////////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename T>
SymmetricMatrix<T>& SymmetricMatrix<T>::operator=(const SymmetricMatrix<T>& other)
{
 if (data.size()!=0)
 {
  for (indextype r=0;r<data.size();r++)
   data[r].clear();
 }
 
 ((JMatrix<T> *)this)->operator=((const JMatrix<T> &)other);
 data.resize(this->nr);
 for (indextype r=0;r<this->nr;r++)
 {
     data[r].resize(r+1);
     copy(other.data[r].begin(),other.data[r].end(),data[r].begin());
 }
 return *this;
}

TEMPLATES_OPERATOR(SymmetricMatrix,=)

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Transpose operator, denoted as != in FullMatrix, is not implemented. It does not make sense to transpose a symmetric matrix.

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Constructor to read from a csv file.
// This is implemented by request of a user, but it is strange to store a symmetric matrix in a .csv file, since it will
// use twice the space it really needs...

template <typename T>
SymmetricMatrix<T>::SymmetricMatrix(std::string fname,unsigned char vtype,char csep) : JMatrix<T>(fname,MTYPESYMMETRIC,vtype,csep)
{
    std::string line;
    // This is just to know number of rows
    this->nr=0;
    while (!this->ifile.eof())
    {
        getline(this->ifile,line);
        if (!this->ifile.eof())
            this->nr++;
    }
    if (this->nr != this->nc)
    {
        std::string err = "csv table in file "+fname+" has different number of rows and columns (as inferred from its header).\n";
        err += "   It is not square, so it cannot be stored as a symmetric matrix.\n";
        Rcpp::stop(err);
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
        Rcpp::Rcout << "WARNING: you are trying to read a symmetric matrix from a .csv file. You .csv file MUST contain a square matrix,\n";
        Rcpp::Rcout << "         but only the lower-triangular matrix (incuding the main diagonal) of it will be stored. Values at the\n";
        Rcpp::Rcout << "         upper-triangular matrix will be read just to check the number of them and immediately ignored.\n";
    }
    
    data.resize(this->nr);
    for (indextype r=0;r<this->nr;r++)
    {
     data[r].resize(r+1);
     data[r].assign(r+1,T(0));
    }
    
    // Reposition pointer at the beginnig and re-read first (header) line
    // and no, seekg does not work (possibly, because we had reached the eof ???)
    this->ifile.close();
    this->ifile.open(fname.c_str());
    size_t p=0;
    getline(this->ifile,line);
    // No need to process first line here, it was done at the parent's class constructor
   
    if (DEB & DEBJM)
        Rcpp::Rcout << "Reading line... ";
    while (!this->ifile.eof())
    {
        if ( (DEB & DEBJM) && !(p%1000))
        {
            Rcpp::Rcout << p << " ";
            Rcpp::Rcout.flush();
        }
        getline(this->ifile,line);
        if (!this->ifile.eof())
        {
          if (!this->ProcessDataLineCsvForSymmetric(line,csep,p,data[p]))
          {
              std::ostringstream errst;
              errst << "Format error reading line " << p << " of file " << fname << ".\n";
              Rcpp::stop(errst.str());
          }
          p++;
          
          if ( (DEB & DEBJM) && (this->nr>1000) && (!(p % 100)) )
           Rcpp::Rcout << p << " "; 
        }
    }
    
    if (DEB & DEBJM)
    {
        Rcpp::Rcout << "\nRead " << p << " data lines of file " << fname;
        if (p==this->nr)
            Rcpp::Rcout << ", as expected.\n";
        else
            Rcpp::Rcout << " instead of " << this->nr << ".\n";
    }
    
    // No call to ReadMetadata must be done here, since these data are NOT binary. The column names were read by the parent's class constructor
    // and the row names are stored as they are read by the former call to ProcessDataLine
    
    
    this->ifile.close();  
}

TEMPLATES_CONST(SymmetricMatrix,SINGLE_ARG(std::string fname,unsigned char vtype,char csep))

//////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename T>
bool SymmetricMatrix<T>::TestDistDisMat()
{
    bool ret=true;
    
    // Check main diagonal. All elements must be 0
    indextype t=0;
    while (ret && t<this->nr)
    {
        ret = (data[t][t]==T(0));
        t++;
    }
    if (!ret)
    {
        Rcpp::Rcerr << "Element (" << t << "," << t << ") and possibly others is/are not 0.\n";
        return false;
    }
    
    // Now, check element outside main diagonal. All elements must be strictly possitive.
    t=1;
    indextype t2;
    while (ret && t<this->nr)
    {
        t2=0;
        while (ret && t2<t)
        {
            ret = (data[t][t2]>=T(0));
            t2++;
        }
        if (!ret)
        {
            t2--;
            Rcpp::Rcerr << "Element (" << t << "," << t2 << ") and possibly others is/are negative, indeed it is " << data[t][t2] << "\n";
            return false;
        }
        t++;
    }
    
    return ret;
}

TEMPLATES_FUNC(bool,SymmetricMatrix,TestDistDisMat,)

//////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef WITH_CHECKS_MATRIX
template <typename T>
T SymmetricMatrix<T>::Get(indextype r,indextype c)
{ 

    if ((r>=this->nr) || (c>=this->nc))
    {
    	std::ostringstream errst;
        errst << "Runtime error in SymmetricMatrix<T>::Get: at least one index (" << r << " or " << c << ") out of bounds.\n";
        errst << "This matrix was of dimension (" << this->nr << " x " << this->nc << ")\n";
        Rcpp::stop(errst.str());
    }
    return (c<=r) ? data[r][c] : data[c][r];
}

TEMPLATES_FUNCR(SymmetricMatrix,Get,SINGLE_ARG(indextype r,indextype c))
#endif

////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef WITH_CHECKS_MATRIX
template <typename T>
void SymmetricMatrix<T>::Set(indextype r,indextype c,T v)
{
    if ((r>=this->nr) || (c>=this->nc))
    {
    	std::ostringstream errst;
        errst << "Runtime error in SymmetricMatrix<T>::Set: at least one index (" << r << " or " << c << ") out of bounds.\n";
        errst << "This matrix was of dimension (" << this->nr << " x " << this->nc << ")\n";
        Rcpp::stop(errst.str());
    }
    if (c<=r) 
        data[r][c]=v;
    else 
        data[c][r]=v;
}

TEMPLATES_SETFUNC(void,SymmetricMatrix,Set,SINGLE_ARG(indextype r,indextype c),v)
#endif

////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename T>
T SymmetricMatrix<T>::GetRowSum(indextype r)
{
#ifdef WITH_CHECKS_MATRIX
    if (r>=this->nr)
    {
    	std::ostringstream errst;
        errst << "Runtime error in SymmetricMatrix<T>::GetRowSum: row index " << r << " out of bounds.\n";
        errst << "This matrix was of dimension (" << this->nr << " x " << this->nr << ")\n";
        Rcpp::stop(errst.str());
    }
#endif
 T sum=T(0);
 for (indextype c=0; c<this->nc; c++)
     sum += ((c<=r) ? data[r][c] : data[c][r]);
 return sum;
}

TEMPLATES_FUNCR(SymmetricMatrix,GetRowSum,indextype r)

//////////////////////////////////////////////////////////////////////////////////////////////////////
template <typename T>
void SymmetricMatrix<T>::WriteBin(std::string fname)
{
    ((JMatrix<T> *)this)->WriteBin(fname,MTYPESYMMETRIC);
    
    if (DEB & DEBJM)
    {
     Rcpp::Rcout << "Writing binary matrix " << fname << std::endl;
     Rcpp::Rcout.flush();
    }
    
    T *ddata = new T [this->nr];
    
    for (indextype r=0;r<this->nr;r++)
    {
        for (indextype c=0;c<=r;c++)
            ddata[c]=data[r][c];
        this->ofile.write((const char *)ddata,(r+1)*sizeof(T));
    }
    
    delete[] ddata;
    
    unsigned long long endofbindata = this->ofile.tellp();
    
    if (DEB & DEBJM)
     Rcpp::Rcout << "End of block of binary data at offset " << endofbindata << "\n";
     
    this->WriteMetadata();                // Here we must write the metadata at the end of the binary contents of the matrix
    
    this->ofile.write((const char *)&endofbindata,sizeof(unsigned long long));  // This writes the point where binary data ends at the end of the file
    
    this->ofile.close();
}

TEMPLATES_FUNC(void,SymmetricMatrix,WriteBin,std::string fname)

/////////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename T>
void SymmetricMatrix<T>::WriteCsv(std::string fname,char csep,bool withquotes)
{
    // Remember: this writes the header, even if there were no column names (then, the header will be "","C1","C2",...)
    ((JMatrix<T> *)this)->WriteCsv(fname,csep,withquotes);

    // Header has been written (unless the matrix had no columns...)
    // But if it would have columns, but no rows, there is nothing after the header. The .csv would have a single line (the header)
    if ((this->nc==0) || (this->nr==0))
    {
     this->ofile.close();
     return;
    }

    int p = std::numeric_limits<T>::max_digits10;

    // We have rows to write; otherwise we would have returned four lines ago...
    indextype rns=this->rownames.size();
    
    for (indextype r=0;r<this->nr;r++)
    {
        if (rns>0)
            this->ofile << FixQuotes(this->rownames[r],withquotes) << csep;
        else
        {
         if (withquotes)
            this->ofile << "\"R" << r+1 << "\"";
         else
            this->ofile << "R" << r+1;
         this->ofile << csep;   // Blank empty field at the beginning of each line
        }

        for (indextype c=0;c<r+1;c++)
            this->ofile << std::setprecision(p) << data[r][c] << csep;
            
        // Csv version writes all elements, even those not physically stored.
        for (indextype c=r+1;c<this->nr-1;c++)
            this->ofile << std::setprecision(p) << data[c][r] << csep;
        this->ofile << std::setprecision(p) << data[this->nr-1][r] << std::endl;
    }
    this->ofile.close();
}

TEMPLATES_FUNC(void,SymmetricMatrix,WriteCsv,SINGLE_ARG(std::string fname,char csep,bool withquotes))

///////////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename T>
float SymmetricMatrix<T>::GetUsedMemoryMB()
{
    unsigned long long num_elem=(unsigned long long)this->nr*(unsigned long long)(this->nr+1)/2;
    Rcpp::Rcout << num_elem << " elements of " << sizeof(T) << " bytes each with accounts for ";
    return float(num_elem)*float(sizeof(T))/(1024.0*1024.0);
}

TEMPLATES_FUNC(float,SymmetricMatrix,GetUsedMemoryMB,)

