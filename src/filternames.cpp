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

/***********************************************************************************
 *
 * The set of functions to filter rows or columns from a jmatrix file.
 *
 ***********************************************************************************/
 
#include  <Rcpp.h>

#include <fullmatrix.h>
#include <sparsematrix.h>
#include <symmetricmatrix.h>
#include <cmath>

extern unsigned char DEB;

// Auxiliary functions to filter genes

std::vector<std::string> FilterAndCheckNames(std::vector<std::string> &names,std::vector<std::string> &gnames,bool byrows,std::vector<bool> &remain,
					     indextype oldnum,indextype &newnrows,indextype &newncols)
{
 std::vector<std::string> remainnames;
 size_t j;
 std::vector<int> timesfound(gnames.size()); 
 for (size_t i=0; i<gnames.size(); i++)
  timesfound[i]=0;
  
 for (size_t i=0; i<names.size(); i++)
 {
  j=0;
  while ( j<gnames.size() && names[i]!=gnames[j] )
   j++;
  if (j<gnames.size())
  {
   timesfound[j]++;
   remain.push_back(true);
   remainnames.push_back(names[i]);
  }
  else
   remain.push_back(false);
 }
 
 // Just for checking... 
 bool anynotfound=false;
 bool anyrepeated=false;
 std::string warn1="";
 std::string warn2="These names were in the passed list of names to be kept and they were more than once in the list of "+(byrows ? std::string("row") : std::string("column"))+" names of the original matrix:\n";
 for (size_t i=0; i< gnames.size(); i++)
 {
  if (timesfound[i]==0)
  { 
   anynotfound = true;
   warn1 += (gnames[i]+" ");
  }
  if (timesfound[i]>1)
  {
   anyrepeated = true;
   warn2 += (gnames[i]+" ");
  } 
 }
 
 if (anynotfound)
 {
  if (byrows)
   Rcpp::warning("These names were in the passed list of names to be kept but they were not in the list of row names of the original matrix:\n");
  else
   Rcpp::warning("These names were in the passed list of names to be kept but they were not in the list of column names of the original matrix:\n");
  Rcpp::Rcerr << warn1 << "\n";
 }
   
 if (anyrepeated)
 {
  if (byrows)
   Rcpp::warning("These names were in the passed list of names to be kept and they were more than once in the list of row names of the original matrix:\n");
  else
   Rcpp::warning("These names were in the passed list of names to be kept and they were more than once in the list of column names of the original matrix:\n");
  Rcpp::Rcerr << warn2 << "\n"; 
 }
 
 if (byrows)
 {
  newnrows=remainnames.size();
  newncols=oldnum;
  if (newnrows==0)
   Rcpp::stop("After filtering, no rows remain in the matrix. Nothing to be saved.\n");
  if (DEB & DEBSC)
   Rcpp::Rcout << newnrows << " rows of the " << names.size() << " in the original matrix will be kept.\n";
 }
 else
 {
  newnrows=oldnum;
  newncols=remainnames.size();
  if (newncols==0)
   Rcpp::stop("After filtering, no columns remain in the matrix. Nothing to be saved.\n");
  if (DEB & DEBSC)
   Rcpp::Rcout << newncols << " columns of the " << names.size() << " in the original matrix will be kept.\n";
 }
 return(remainnames);
}

template <typename T>
void FilterF(FullMatrix<T> &M,std::vector<std::string> gnames,bool byrows,std::string filname)
{
 std::vector<std::string> names=(byrows ? M.GetRowNames() : M.GetColNames());
 indextype on=(byrows ? M.GetNCols() : M.GetNRows());
 std::vector<bool> remain;
 indextype newnrows,newncols;
 std::vector<std::string> remainnames=FilterAndCheckNames(names,gnames,byrows,remain,on,newnrows,newncols); 
 FullMatrix<T> Rem(newnrows,newncols);
 if (byrows)
 {
  indextype nr=0;
  for (indextype r=0;r<M.GetNRows();r++)
   if (remain[r])
   {
    for (indextype c=0;c<M.GetNCols();c++)
     Rem.Set(nr,c,M.Get(r,c));
    nr++;
   }
  Rem.SetRowNames(remainnames);
  Rem.SetColNames(M.GetColNames());
 }
 else
 {
  indextype nc=0;
  for (indextype c=0;c<M.GetNCols();c++)
   if (remain[c])
   {
    for (indextype r=0;r<M.GetNRows();r++)
     Rem.Set(r,nc,M.Get(r,c));
    nc++;
   }
  Rem.SetRowNames(M.GetRowNames());
  Rem.SetColNames(remainnames);
 }
 Rem.SetComment(M.GetComment());
 Rem.WriteBin(filname);
}

template <typename T>
void FilterS(SparseMatrix<T> &M,std::vector<std::string> gnames,bool byrows,std::string filname)
{
 std::vector<std::string> names=(byrows ? M.GetRowNames() : M.GetColNames());
 indextype on=(byrows ? M.GetNCols() : M.GetNRows());
 std::vector<bool> remain;
 indextype newnrows,newncols;
 std::vector<std::string> remainnames=FilterAndCheckNames(names,gnames,byrows,remain,on,newnrows,newncols);

 SparseMatrix<T> Rem(newnrows,newncols);
 if (byrows)
 {
  indextype nr=0;
  for (indextype r=0;r<M.GetNRows();r++)
   if (remain[r])
   {
    for (indextype c=0;c<M.GetNCols();c++)
     Rem.Set(nr,c,M.Get(r,c));
    nr++;
   }
  Rem.SetRowNames(remainnames);
  Rem.SetColNames(M.GetColNames());
 }
 else
 {
  indextype nc=0;
  for (indextype c=0;c<M.GetNCols();c++)
   if (remain[c])
   {
    for (indextype r=0;r<M.GetNRows();r++)
     Rem.Set(r,nc,M.Get(r,c));
    nc++;
   }
  Rem.SetRowNames(M.GetRowNames());
  Rem.SetColNames(remainnames);
 }
 Rem.SetComment(M.GetComment());
 Rem.WriteBin(filname);
}

void FilterAndSaveFull(std::string fname,unsigned char ctype,bool genesbyrows,std::vector<std::string> gnames,std::string filname)
{
 switch (ctype)
    {
        case UCTYPE: { FullMatrix<unsigned char> M(fname); FilterF<unsigned char>(M,gnames,genesbyrows,filname); } break;
        case SCTYPE: { FullMatrix<char> M(fname); FilterF<char>(M,gnames,genesbyrows,filname); } break; 
        case USTYPE: { FullMatrix<unsigned short> M(fname); FilterF<unsigned short>(M,gnames,genesbyrows,filname); } break; 
        case SSTYPE: { FullMatrix<short> M(fname); FilterF<short>(M,gnames,genesbyrows,filname); } break; 
        case UITYPE: { FullMatrix<unsigned int> M(fname); FilterF<unsigned int>(M,gnames,genesbyrows,filname); } break; 
        case SITYPE: { FullMatrix<int> M(fname); FilterF<int>(M,gnames,genesbyrows,filname); } break; 
        case ULTYPE: { FullMatrix<unsigned long> M(fname); FilterF<unsigned long>(M,gnames,genesbyrows,filname); } break; 
        case SLTYPE: { FullMatrix<long> M(fname); FilterF<long>(M,gnames,genesbyrows,filname); } break; 
        case FTYPE:  { FullMatrix<float> M(fname); FilterF<float>(M,gnames,genesbyrows,filname); } break;  
        case DTYPE:  { FullMatrix<double> M(fname); FilterF<double>(M,gnames,genesbyrows,filname); } break;  
        case LDTYPE: { FullMatrix<long double> M(fname); FilterF<long double>(M,gnames,genesbyrows,filname); } break;
        default: Rcpp::stop("Matrix in input file is on unknown data type. Was it created by package jmatrix/parallelpam/scellpam?\n"); break;
    }
}

void FilterAndSaveSparse(std::string fname,unsigned char ctype,bool genesbyrows,std::vector<std::string> gnames,std::string filname)
{
 switch (ctype)
    {
        case UCTYPE: { SparseMatrix<unsigned char> M(fname); FilterS<unsigned char>(M,gnames,genesbyrows,filname); } break;
        case SCTYPE: { SparseMatrix<char> M(fname); FilterS<char>(M,gnames,genesbyrows,filname); } break; 
        case USTYPE: { SparseMatrix<unsigned short> M(fname); FilterS<unsigned short>(M,gnames,genesbyrows,filname); } break; 
        case SSTYPE: { SparseMatrix<short> M(fname); FilterS<short>(M,gnames,genesbyrows,filname); } break; 
        case UITYPE: { SparseMatrix<unsigned int> M(fname); FilterS<unsigned int>(M,gnames,genesbyrows,filname); } break; 
        case SITYPE: { SparseMatrix<int> M(fname); FilterS<int>(M,gnames,genesbyrows,filname); } break; 
        case ULTYPE: { SparseMatrix<unsigned long> M(fname); FilterS<unsigned long>(M,gnames,genesbyrows,filname); } break; 
        case SLTYPE: { SparseMatrix<long> M(fname); FilterS<long>(M,gnames,genesbyrows,filname); } break; 
        case FTYPE:  { SparseMatrix<float> M(fname); FilterS<float>(M,gnames,genesbyrows,filname); } break;  
        case DTYPE:  { SparseMatrix<double> M(fname); FilterS<double>(M,gnames,genesbyrows,filname); } break;  
        case LDTYPE: { SparseMatrix<long double> M(fname); FilterS<long double>(M,gnames,genesbyrows,filname); } break;
        default: Rcpp::stop("Matrix in input file is on unknown data type. Was it created by package jmatrix/parallelpam/scellpam?\n"); break;
    }
}

//' FilterJMatByName
//'
//' Takes a jmatrix binary file containing a table with rows and columns and filters it by name, eliminating the rows or columns whose whose names are not in certain list
//'
//' If the table has no list of names in the requested dimension (rows or colums), an error is rised.\cr
//' The row or column names whose names are not found obviosuly cannot remain, and the program rises a warning indicating for which row/column names this happens.\cr
//' The matrix contained in the filtered file will have the same nature (full or sparse) and the same data type as the original.\cr
//' This function can be used to filter either by row or by column name, with appropriate usage of parameter namesat
//'
//' @param fname   A string with the file name of the original table
//' @param Gn      A list of R strings with the names of the rows or columns that must remain. All others will be filtered out
//' @param filname A string with the file name of the filtered table
//' @param namesat The string "rows" or "cols" indicating if the searched names are in the rows or in the columns of the original table. Default: "rows"
//' @return        No return value, called for side effects (creates a file)
//' @examples
//' Rf <- matrix(runif(48),nrow=6)
//' rownames(Rf) <- c("A","B","C","D","E","F")
//' colnames(Rf) <- c("a","b","c","d","e","f","g","h")
//' tmpfile1=paste0(tempdir(),"/Rfullfloat.bin")
//' tmpfile2=paste0(tempdir(),"/Rfullfloatrowfilt.bin")
//' tmpfile3=paste0(tempdir(),"/Rfullfloatrowcolfilt.bin")
//' tmpcsvfile1=paste0(tempdir(),"/Rfullfloat.csv")
//' tmpcsvfile3=paste0(tempdir(),"/Rfullfloatrowcolfilt.csv")
//' JWriteBin(Rf,tmpfile1,dtype="float",dmtype="full",comment="Full matrix of floats")
//' # Let's keep only rows A, C and E
//' FilterJMatByName(tmpfile1,c("A","C","E"),tmpfile2,namesat="rows")
//' # and from the result, let's keep only columns b, d and g
//' FilterJMatByName(tmpfile2,c("b","d","g"),tmpfile3,namesat="cols")
//' JMatToCsv(tmpfile1,tmpcsvfile1)
//' JMatToCsv(tmpfile3,tmpcsvfile3)
//' # You can now compare both ASCII/csv files
//' @export
// [[Rcpp::export]]
void FilterJMatByName(std::string fname,Rcpp::StringVector Gn,std::string filname,std::string namesat="rows")
{
 if ( (namesat != "rows") && (namesat != "cols") && (namesat != "columns") )
  Rcpp::stop("Valid values for parameter namesat are only 'rows' and 'cols'.\n");
  
 unsigned char mtype,ctype,endian,mdinfo;
 indextype nrows,ncols;
 
 MatrixType(fname,mtype,ctype,endian,mdinfo,nrows,ncols);
 
 std::vector<std::string> names;
 for (long int i=0; i< Gn.length(); i++)
  names.push_back(std::string(Gn[i]));
  
 switch (mtype)
 {
        case MTYPEFULL:         FilterAndSaveFull(fname,ctype,(namesat=="rows"),names,filname); break;
        case MTYPESPARSE:       FilterAndSaveSparse(fname,ctype,(namesat=="rows"),names,filname); break;
        case MTYPESYMMETRIC:    Rcpp::stop("This function cannot be applied to symmetric matrices, only to full or sparse matrices.\n"); break;
        default:                Rcpp::stop("Unknown matrix type. Was the input file generated by the jmatrix/parallelpam/scellpam packages?\n"); break;
 }
}

