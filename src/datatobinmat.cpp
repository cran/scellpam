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
 * The set of functions to import/export data form/to the most used formats:
 * CSV (comma separated values), Seurat/dgCMatrix and Sce (Single cell experiment).
 *
 * The first one does not need anything else, since it is self-described.
 *
 * The Sce might be helped by the use of package BiocGenerics, concretely of its
 * function 'counts' which extracts the matrix of counts because it knows the internal
 * structure of some other objects in several packages that use it. It will work, as
 * long as BiocGenerics follows timely the changes in object structure that the designers
 * of the particular packages may introduce in the future. This is why we have opted
 * by supporting Sce asking the user to extract the data matrix from the Sce object.
 *
 * Similarly, the Seurat format cannot be treated globally because the internal
 * structure of the objects is not unified, so again the counts must be extracted
 * from within the object. Fortunately, it seems that all objects compatible with
 * Seurat store their counts as a dgCMatrix (sparse matrix). That is why we provide
 * a function that gets this type of matrices. 
 *
 ***********************************************************************************/
 
#include  <Rcpp.h>

#include <fullmatrix.h>
#include <sparsematrix.h>
#include <symmetricmatrix.h>
#include <cmath>

extern unsigned char DEB;

/***********************************************************
 * 
 * Functions related with Csv format import/export
 *
 **********************************************************/
 
// Underlying internal functions


template <typename T>
void CsvDataToBinMat(std::string ifname,std::string ofname,unsigned char vtype,std::string ctype,char csep,bool isfull,bool transpose,std::string comment)
{
 if (isfull)
 {
  FullMatrix<T> M(ifname,vtype,csep);
  
  if (ctype!="raw")
   M.SelfColNorm(ctype);      // In this case, the .csv must be read as is; therefore, the normalization is done as this functions says it does: by columns.
  
  if (comment!="")
   M.SetComment(comment);
 
  if (!transpose)
   M.WriteBin(ofname);
  else
  {
   FullMatrix<T> Mt;
 
   // I had to use this as transpose-assignment operator. Postfix operators (to use something like Mt = M') are not allowed in C++
   Mt != M;
  
   Mt.WriteBin(ofname);
  }
 }
 else
 {
  SparseMatrix<T> M(ifname,vtype,csep);
 
  if (ctype!="raw")
   M.SelfColNorm(ctype);      // In this case, the .csv must be read as is; therefore, the normalization is done as this functions says it does: by columns.
   
  if (comment!="")
   M.SetComment(comment);
 
  if (!transpose)
   M.WriteBin(ofname);
  else
  {
   SparseMatrix<T> Mt;
  
   // I had to use this as transpose-assignment operator. Postfix operators (to use something like Mt = M') are not allowed in C++
   Mt != M;
  
   Mt.WriteBin(ofname);
  }
 }
}

//' CsvToJMat
//'
//' Gets a csv/tsv file and writes to a disk file the binary matrix of counts contained in it in the jmatrix binary format.\cr
//' First line of the .csv is supposed to have the field names.\cr
//' First column of each line is supposed to have the field name.\cr
//' The fields are supposed to be separated by one occurrence of a character-field sepparator (usually, comma or tab)
//' .tsv files can be read with this function, too, setting the csep argument to '\\t'
//'
//' The parameter transpose has the default value of FALSE. But don't forget to set it to TRUE if you want the cells
//' (which in single cell common practice are by columns) to be written by rows. This will be needed later to calculate
//' the dissimilarity matrix, if this is the next step of your workflow. See help of CalcAndWriteDissimilarityMatrix
//' 
//' @param ifname    A string with the name of the .csv/.tsv text file.
//' @param ofname    A string with the name of the binary output file.
//' @param mtype     A string to indicate the matrix type: 'full' or 'sparse'. Default: 'sparse'
//' @param csep      The character used as separator in the .csv file. Default: ',' (comma) (Set to '\\t' for .tsv)
//' @param ctype     The string 'raw' or 'log1' to write raw counts or log(counts+1), or the normalized versions, 'rawn' and 'log1n', which normalize ALWAYS BY COLUMNS (before transposition, if requested to transpose). The logarithm is taken base 2. Default: raw
//' @param valuetype The data type to store the matrix. It must be one of the strings 'uint32', 'float' or 'double'. Default: float
//' @param transpose Boolean to indicate if the matrix should be transposed before writing. See Details for a comment about this. Default: FALSE
//' @param comment   A comment to be stored with the matrix. Default: "" (no comment)
//' @return          No return value, called for side effects (creates a file)
//' @examples
//' # Since we have no a .csv file to test, we will generate one with another funcion of this package
//' Rf <- matrix(runif(48),nrow=6)
//' rownames(Rf) <- c("A","B","C","D","E","F")
//' colnames(Rf) <- c("a","b","c","d","e","f","g","h")
//' tmpfile1=paste0(tempdir(),"/Rfullfloat.bin")
//' tmpfile2=paste0(tempdir(),"/Rfullfloat2.bin")
//' tmpcsvfile1=paste0(tempdir(),"/Rfullfloat.csv")
//' JWriteBin(Rf,tmpfile1,dtype="float",dmtype="full",comment="Full matrix of floats")
//' JMatToCsv(tmpfile1,tmpcsvfile1)
//' CsvToJMat(tmpcsvfile1,tmpfile2)
//' # It can be checked that files Rfullfloat.bin and Rfullfloat2.bin contain the same data
//' # (even they differ in the comment, which has been eliminated when converting to csv)
//' @export
// [[Rcpp::export]]
void CsvToJMat(std::string ifname, std::string ofname,std::string mtype = "sparse",char csep=',', std::string ctype = "raw", std::string valuetype = "float", bool transpose = false,std::string comment = "")
{
 if ((ctype != "raw") && (ctype != "log1") && (ctype != "rawn") && (ctype != "log1n"))
 {
  Rcpp::stop("The ctype argument can take only one of the string values 'raw', 'log1', 'rawn' and 'log1n'\n");
  return;
 }
 if ((mtype != "full") && (mtype != "sparse"))
 {
  Rcpp::stop("The mtype argument can take only one of the string values 'full' or 'sparse'\n");
  return;
 }
 bool isfull=(mtype=="full");
 if ((valuetype != "float") && (valuetype != "double") && (valuetype != "uint32"))
 {
  Rcpp::stop("The valuetype argument can take only one of the string values 'uint32', 'float' or 'double'\n");
  return;
 }
 if ((valuetype == "uint32") && (ctype == "log1"))
 {
  Rcpp::stop("Rescaling as log(counts+1) requires output type to be float or double, not uint32.\n");
  return;
 }
 
 if (valuetype=="uint32")
  CsvDataToBinMat<unsigned int>(ifname,ofname,ULTYPE,ctype,csep,isfull,transpose,comment);
  
 if (valuetype=="float")
  CsvDataToBinMat<float>(ifname,ofname,FTYPE,ctype,csep,isfull,transpose,comment);
  
 if (valuetype=="double")
  CsvDataToBinMat<double>(ifname,ofname,DTYPE,ctype,csep,isfull,transpose,comment);
}

//' JMatToCsv
//'
//' Writes a binary matrix in the jmatrix package format as a .csv file. This is mainly for checking/inspection and
//' to load the data from R as read.table, if the memory of having all data as doubles allows doing such thing.
//'
//' @param ifile      String with the file name that contains the binary data.
//' @param csvfile    String with the file name that will contain the data as csv.
//' @param csep       Character used as separator. Default: , (comma)
//' @param withquotes boolean to mark if row and column names in the .csv file must be written surrounded by doble quotes. Default: FALSE
//' @return           No return value, called for side effects (creates a file)
//' @examples
//' Rf <- matrix(runif(48),nrow=6)
//' rownames(Rf) <- c("A","B","C","D","E","F")
//' colnames(Rf) <- c("a","b","c","d","e","f","g","h")
//' tmpfile1=paste0(tempdir(),"/Rfullfloat.bin")
//' tmpcsvfile1=paste0(tempdir(),"/Rfullfloat.csv")
//' JWriteBin(Rf,tmpfile1,dtype="float",dmtype="full",comment="Full matrix of floats")
//' JMatToCsv(tmpfile1,tmpcsvfile1)
//' @export
// [[Rcpp::export]]
void JMatToCsv(std::string ifile, std::string csvfile, char csep =',',bool withquotes=false)
{
 unsigned char mtype,ctype,endian,mdinf;
 indextype nrows,ncols;
 
 MatrixType(ifile,mtype,ctype,endian,mdinf,nrows,ncols);
 
 if (mtype==MTYPEFULL)
 {
     switch (ctype)
     {
        case UCTYPE: { FullMatrix<unsigned char> M(ifile); M.WriteCsv(csvfile,csep,withquotes); break; };
        case SCTYPE: { FullMatrix<char> M(ifile); M.WriteCsv(csvfile,csep,withquotes); break; };
        case USTYPE: { FullMatrix<unsigned short> M(ifile); M.WriteCsv(csvfile,csep,withquotes); break; };
        case SSTYPE: { FullMatrix<short> M(ifile); M.WriteCsv(csvfile,csep,withquotes); break; };
        case UITYPE: { FullMatrix<unsigned int> M(ifile); M.WriteCsv(csvfile,csep,withquotes); break; };
        case SITYPE: { FullMatrix<int> M(ifile); M.WriteCsv(csvfile,csep,withquotes); break; };
        case ULTYPE: { FullMatrix<unsigned long> M(ifile); M.WriteCsv(csvfile,csep,withquotes); break; };
        case SLTYPE: { FullMatrix<long> M(ifile); M.WriteCsv(csvfile,csep,withquotes); break; };
        case FTYPE:  { FullMatrix<float> M(ifile); M.WriteCsv(csvfile,csep,withquotes); break; };
        case DTYPE:  { FullMatrix<double> M(ifile); M.WriteCsv(csvfile,csep,withquotes); break; };
        case LDTYPE: { FullMatrix<long double> M(ifile); M.WriteCsv(csvfile,csep,withquotes); break; };
        default: break;
    }
 }
 if (mtype==MTYPESPARSE)
 {
     switch (ctype)
     {
        case UCTYPE: { SparseMatrix<unsigned char> M(ifile); M.WriteCsv(csvfile,csep,withquotes); break; };
        case SCTYPE: { SparseMatrix<char> M(ifile); M.WriteCsv(csvfile,csep,withquotes); break; };
        case USTYPE: { SparseMatrix<unsigned short> M(ifile); M.WriteCsv(csvfile,csep,withquotes); break; };
        case SSTYPE: { SparseMatrix<short> M(ifile); M.WriteCsv(csvfile,csep,withquotes); break; };
        case UITYPE: { SparseMatrix<unsigned int> M(ifile); M.WriteCsv(csvfile,csep,withquotes); break; };
        case SITYPE: { SparseMatrix<int> M(ifile); M.WriteCsv(csvfile,csep,withquotes); break; };
        case ULTYPE: { SparseMatrix<unsigned long> M(ifile); M.WriteCsv(csvfile,csep,withquotes); break; };
        case SLTYPE: { SparseMatrix<long> M(ifile); M.WriteCsv(csvfile,csep,withquotes); break; };
        case FTYPE:  { SparseMatrix<float> M(ifile); M.WriteCsv(csvfile,csep,withquotes); break; }; 
        case DTYPE:  { SparseMatrix<double> M(ifile); M.WriteCsv(csvfile,csep,withquotes); break; };
        case LDTYPE: { SparseMatrix<long double> M(ifile); M.WriteCsv(csvfile,csep,withquotes); break; };
        default: break;
    }
 }
 if (mtype==MTYPESYMMETRIC)
 {
     switch (ctype)
     {
        case UCTYPE: { SymmetricMatrix<unsigned char> M(ifile); M.WriteCsv(csvfile,csep,withquotes); break; };
        case SCTYPE: { SymmetricMatrix<char> M(ifile); M.WriteCsv(csvfile,csep,withquotes); break; };
        case USTYPE: { SymmetricMatrix<unsigned short> M(ifile); M.WriteCsv(csvfile,csep,withquotes); break; };
        case SSTYPE: { SymmetricMatrix<short> M(ifile); M.WriteCsv(csvfile,csep,withquotes); break; };
        case UITYPE: { SymmetricMatrix<unsigned int> M(ifile); M.WriteCsv(csvfile,csep,withquotes); break; };
        case SITYPE: { SymmetricMatrix<int> M(ifile); M.WriteCsv(csvfile,csep,withquotes); break; };
        case ULTYPE: { SymmetricMatrix<unsigned long> M(ifile); M.WriteCsv(csvfile,csep,withquotes); break; };
        case SLTYPE: { SymmetricMatrix<long> M(ifile); M.WriteCsv(csvfile,csep,withquotes); break; };
        case FTYPE:  { SymmetricMatrix<float> M(ifile); M.WriteCsv(csvfile,csep,withquotes); break; }; 
        case DTYPE:  { SymmetricMatrix<double> M(ifile); M.WriteCsv(csvfile,csep,withquotes); break; };
        case LDTYPE: { SymmetricMatrix<long double> M(ifile); M.WriteCsv(csvfile,csep,withquotes); break; };
        default: break;
    }
 }
}

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
 std::string warn1="These genes were in the passed list of genes to be kept but they were not in the list of "+(byrows ? std::string("row") : std::string("column"))+" names of the original matrix:\n";
 std::string warn2="These genes were in the passed list of genes to be kept and they were more than once in the list of "+(byrows ? std::string("row") : std::string("column"))+" names of the original matrix:\n";
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
  Rcpp::warning(warn1.c_str());
 if (anyrepeated)
  Rcpp::warning(warn2.c_str());
 
 if (byrows)
 {
  newnrows=remainnames.size();
  newncols=oldnum;
  if (DEB & DEBSC)
   Rcpp::Rcout << newnrows << " rows of the " << names.size() << " in the original matrix will be kept.\n";
 }
 else
 {
  newnrows=oldnum;
  newncols=remainnames.size();
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
        default: Rcpp::stop("Matrix in input file is on unknown data type. Was it created by package scellpam?\n"); break;
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
        default: Rcpp::stop("Matrix in input file is on unknown data type. Was it created by package scellpam?\n"); break;
    }
}

//' FilterJMatByName
//'
//' Takes a jmatrix binary file containing a table with cells and genes and filter the genes or cell by name, eliminating those whose names are not in certain list
//'
//' If the table has no list of names in the requested dimension (rows or colums), an error is rised.\cr
//' The gene or cell names whose names are not found obviosuly cannot remain, and the program rises a warning indicating for which gene/cell names this happens.\cr
//' The matrix contained in the filtered file will have the same nature (full or sparse) and the same data type as the original.\cr
//' This function can be used to filter either by cell or by gene name, with appropriate usage of parameter namesat
//'
//' @param fname   A string with the file name of the original table
//' @param Gn      A list of R strings with the names of the genes or cells that must remain. All others will be filtered out
//' @param filname A string with the file name of the filtered table
//' @param namesat The string "rows" or "cols" indicating if the genes/cells at the original table are by rows or by colums. Default: "rows"
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
        default:                Rcpp::stop("Unknown matrix type. Was the input file generated by the scellpam package?\n"); break;
 }
}


/***********************************************************
 * 
 * Functions related with Seurat/dgCMatrix format import/export
 *
 **********************************************************/

// Underlying internal functions

template <typename T>
void PrepareFull(FullMatrix<T> &M,std::string ctype,bool transpose,
                                  Rcpp::StringVector rownames,Rcpp::StringVector colnames,
                                  std::string comment)
{
 if (ctype!="raw")
 {
  if (transpose)
   M.SelfRowNorm(ctype);      // Normalization is done by rows since this matrix is really transposed from it origin
  else
   M.SelfColNorm(ctype);      // Normalization is done as this functions says it does: by columns.
 }
   
 if (comment!="")
  M.SetComment(comment); 
    
 if (DEB & DEBSC)
  Rcpp::Rcout << "Attaching vector of " << colnames.size() << " as " << (transpose ? "row" : "column") << " names and vector of " << rownames.size() << " as " << (transpose ? "row" : "column") << "names.\n";
    
 if (transpose)
 {
  M.SetColNames(rownames);   
  M.SetRowNames(colnames);
 }
 else
 {
  M.SetRowNames(rownames);  
  M.SetColNames(colnames);
 }
}

template <typename T>
void PrepareSparse(SparseMatrix<T> &M,std::string ctype,bool transpose,
                                  Rcpp::StringVector rownames,Rcpp::StringVector colnames,
                                  std::string comment)
{
 if (ctype!="raw")
 {
  if (transpose)
   M.SelfRowNorm(ctype);      // Normalization is done by rows since this matrix is really transposed from it origin
  else
   M.SelfColNorm(ctype);      // Normalization is done as this functions says it does: by columns.
 }
   
 if (comment!="")
  M.SetComment(comment);
 
 if (DEB & DEBSC)
  Rcpp::Rcout << "Attaching vector of " << colnames.size() << " as " << (transpose ? "row" : "column") << " names and vector of " << rownames.size() << " as " << (transpose ? "row" : "column") << "names.\n";
  
 if (transpose)
 {
  M.SetColNames(rownames);   
  M.SetRowNames(colnames);
 }
 else
 {
  M.SetRowNames(rownames);  
  M.SetColNames(colnames);
 }
}

template <typename T>
void dgCMatrixDataToBinMat(std::string fname,std::string ctype,indextype nrows,indextype ncols,bool isfull,bool transpose,
                              Rcpp::NumericVector rowindexes, Rcpp::NumericVector colacc, Rcpp::NumericVector values,
                              Rcpp::StringVector rownames,Rcpp::StringVector colnames,std::string comment)
{
 if (DEB & DEBSC)
 {
  Rcpp::Rcout << "Reading data to put in ";
  Rcpp::Rcout << (transpose ? "transposed" : "non-transposed");
  Rcpp::Rcout << (isfull ? " full" : " sparse");
  Rcpp::Rcout << " matrix. This may be slow. Please, wait...\n";
 } 
  
 // Number of element in current column
 indextype numcpcol;
 // Main index to run through values
 indextype i;
 
 if (isfull)
 {   
  FullMatrix<T> M(ncols,nrows);
        
  i=0; 
  for (indextype col=0;col<ncols;col++)
  {
   numcpcol=colacc[col+1]-colacc[col];
   
   if (transpose)
    for (indextype wrow=0;wrow<numcpcol;wrow++)
     M.Set(col,rowindexes[i+wrow],T(values[i+wrow]));
   else
    for (indextype wrow=0;wrow<numcpcol;wrow++)
     M.Set(rowindexes[i+wrow],col,T(values[i+wrow]));
   
   i+=numcpcol;
  }
   
  PrepareFull(M,ctype,transpose,rownames,colnames,comment);
  M.WriteBin(fname);
 }
 else
 {
  SparseMatrix<T> M(ncols,nrows);
    
  i=0; 
  for (indextype col=0;col<ncols;col++)
  {
   numcpcol=colacc[col+1]-colacc[col];
   
   if (transpose)
    for (indextype wrow=0;wrow<numcpcol;wrow++)
     M.Set(col,rowindexes[i+wrow],T(values[i+wrow]));
   else
    for (indextype wrow=0;wrow<numcpcol;wrow++)
     M.Set(rowindexes[i+wrow],col,T(values[i+wrow]));
     
   i+=numcpcol;
  }
   
  PrepareSparse(M,ctype,transpose,rownames,colnames,comment);
  M.WriteBin(fname);
 }
}

//' dgCMatToJMat
//'
//' Gets a dgCMatrix object (sparse matrix of the `Matrix` package) and writes to a disk file the binary matrix of counts contained in it in the jmatrix binary format. Plase, see Details
//' below to know more about the extraction of the sparse matrices from Seurat or similar single cell formats.
//' 
//' We have found that, in some Seurat objects, the dgCMatrix to be passed to this function can be extracted as q@assays$RNA@counts, being q the Seurat S4 object.\cr
//' In other cases this matrix is obtained as q@raw.data.\cr
//' In any case, we assume that this matrix has slots Dimnames (with a list of strings in Dimnames[[0]] as rownames and Dimnames[[1]]
//' as column names) as long as slots with names i, p and x as described in the documentation of the `Matrix` package on sparse matrices.
//'
//' The parameter transpose has the default value of FALSE. But don't forget to set it to TRUE if you want the cells
//' (which in single cell common practice are by columns) to be written by rows. This will be needed later to calculate
//' the dissimilarity matrix, if this is the next step of your workflow. See help of CalcAndWriteDissimilarityMatrix
//' 
//' 
//' @param q         The dgCMatrix object
//' @param fname     A string with the name of the binary output file
//' @param mtype     A string to indicate the matrix type: 'full' or 'sparse'. Default: 'sparse'
//' @param ctype     The string 'raw' or 'log1' to write raw counts or log(counts+1), or the normalized versions, 'rawn' and 'log1n', which normalize ALWAYS BY COLUMNS (before transposition, if requested to transpose). Default: raw
//' @param valuetype The data type to store the matrix. It must be one of the strings 'uint32', 'float' or 'double'. Default: float
//' @param transpose Boolean to indicate if the matrix should be transposed before writing. See Details for a comment about this. Default: FALSE
//' @param comment   A comment to be stored with the matrix. Default: "" (no comment)
//' @return          No return value, called for side effects (creates a file)
//' @examples
//' # Sorry, we cannot provide an example here, since it would need the load of the Seurat package.
//' # Please, see the vignette for examples
//' @export
// [[Rcpp::export]]
void dgCMatToJMat(Rcpp::S4 q, std::string fname, std::string mtype = "sparse", std::string ctype = "raw", std::string valuetype = "float", bool transpose = false, std::string comment="")
{
 if ((ctype != "raw") && (ctype != "log1") && (ctype != "rawn") && (ctype != "log1n"))
 {
  Rcpp::stop("The ctype argument can take only one of the string values 'raw', 'log1', 'rawn' and 'log1n'\n");
  return;
 }
 if ((mtype != "full") && (mtype != "sparse"))
 {
  Rcpp::stop("The mtype argument can take only one of the string values 'full' or 'sparse'\n");
  return;
 }
 bool isfull=(mtype=="full");
 if ((valuetype != "float") && (valuetype != "double") && (valuetype != "uint32"))
 {
  Rcpp::stop("The valuetype argument can take only one of the string values 'uint32', 'float' or 'double'\n");
  return;
 }
 if ( (valuetype == "uint32") && ((ctype == "log1") || (ctype == "log1n")) )
 {
  Rcpp::stop("Rescaling as log(counts+1) requires output type to be float or double, not uint32, since it produces decimal numbers.\n");
  return;
 }
 
 if ( (valuetype == "uint32") && (ctype == "rawn") )
 {
  Rcpp::stop("Rescaling raw data by normalization requires output type to be float or double, not uint32, since it produces decimal numbers.\n");
  return;
 }
 
 Rcpp::NumericVector dim = q.slot("Dim");
 unsigned int nrows=dim[0];
 unsigned int ncols=dim[1];
 
 Rcpp::List q4 = q.slot("Dimnames");
 Rcpp::StringVector rownames=q4[0];
 Rcpp::StringVector colnames=q4[1];

 if (nrows!=(unsigned int)rownames.length())
  Rcpp::stop("Strange dgCMatrix object. The number of rows in the matrix differs from the length of the vector of row names.\n");
 if (ncols!=(unsigned int)colnames.length())
  Rcpp::stop("Strange dgCMatrix object. The number of columns in the matrix differs from the length of the vector of column names.\n"); 
  
 if (DEB & DEBSC)
 {
  Rcpp::Rcout << "The matrix of counts is of dimension [" << nrows << " x " << ncols << "].\n";
  Rcpp::Rcout << "Row names from " << rownames[0] << " to " << rownames[nrows-1] << "\n";
  Rcpp::Rcout << "Column names from " << colnames[0] << " to " << colnames[ncols-1] << "\n";
 }
  
 if (nrows<2 || ncols<2)
 {
  Rcpp::stop("Count matrix must have at least two rows and two columns.\n");
  return;
 }
 Rcpp::NumericVector rowindexes = q.slot("i");
 Rcpp::NumericVector colacc = q.slot("p");
 Rcpp::NumericVector values = q.slot("x");
 
 if (values.length() < 2) 
 {
  Rcpp::stop("The vector of values in sparse matrix is too small (length is 0 or 1)\n");
  return;
 }
 
 if ((unsigned long long)values.length() > (unsigned long long)nrows*(unsigned long long)ncols)
 {
  Rcpp::stop("The vector of values in sparse matrix has incorrect length (it is too big, bigger than the product of number of rows x number of columns().\n");
  return;
 }
 
 if (rowindexes.length() != values.length())
 {
  Rcpp::stop("The vectors of row indexes and of values in the sparse matrix have not the same length.\n");
  return;
 }
 
 if ((unsigned int)colacc.length() != ncols+1)
 {
  Rcpp::stop("The vector of column accumulated indexes in sparse matrix has not the number of columns plus one as its length.\n");
  return;
 }
 
 if (valuetype=="uint32")
  dgCMatrixDataToBinMat<unsigned int>(fname,ctype,nrows,ncols,isfull,transpose,rowindexes,colacc,values,rownames,colnames,comment);
 
 if (valuetype=="float")
  dgCMatrixDataToBinMat<float>(fname,ctype,nrows,ncols,isfull,transpose,rowindexes,colacc,values,rownames,colnames,comment);
  
 if (valuetype=="double") 
  dgCMatrixDataToBinMat<double>(fname,ctype,nrows,ncols,isfull,transpose,rowindexes,colacc,values,rownames,colnames,comment);
}

//' GetSeuratGroups
//'
//' Returns a numeric vector of integers with the numeric identifier of the group to which each cell in a Seurat object belongs to, if the cells come from 
//' different groups/samples. These numeric identifiers go from 1 to the number of groups; names of original factors are not kept.
//'
//' If q is the Seurat object, this function assumes that\cr
//'  \cr
//' q@meta.data$orig.ident\cr
//'  \cr
//' is the integer vector with this information. We don't know if this is assumed by all software which uses Seurat (probably, not)
//' so this function is likely NOT to work in most cases and therefore is provided just as a convenience that can generate the parameter gr
//' for the BuildAbundanceMatrix. But if the data you have got does not follow these conventions, please don't blame us...
//'
//' @param q		The S4 Seurat object (for example, returned by a call to readRDS('file.rds') where the rds file was written by Seurat).
//' @return     	The numeric integer vector with as many components as cells.
//' @examples
//' # Sorry, we cannot provide an example here, since it would need the load of the Seurat package.
//' # Please, see the vignette for examples
//' @export
// [[Rcpp::export]]
Rcpp::IntegerVector GetSeuratGroups(Rcpp::S4 q)
{
 Rcpp::DataFrame q1 = static_cast<Rcpp::DataFrame>(q.slot("meta.data")); 
 Rcpp::IntegerVector q2 = q1["orig.ident"];
 Rcpp::IntegerVector ret(q2.length());
 
 int min=q2.length();
 int max=0;
 for (indextype i=0; i<(unsigned int)q2.length(); i++)
 {
  if (q2[i]>max)
   max=q2[i];
  if (q2[i]<min)
   min=q2[i];
  ret[i]=q2[i];
 }
 if ((min!=1) || (max==min))
  Rcpp::stop("The minimum group number is not 1, or the maximum group number is 1.\n");
  
 Rcpp::List qt1 = q.slot("assays");
 Rcpp::S4 qt2 = qt1["RNA"];
 Rcpp::S4 qt3 = qt2.slot("counts");
  
 Rcpp::NumericVector dim = qt3.slot("Dim");
 
 Rcpp::List qt4 = qt3.slot("Dimnames");
 Rcpp::StringVector rownames=qt4[0];
 Rcpp::StringVector colnames=qt4[1];

 if (colnames.length()!=q2.length())
  Rcpp::stop("The lengths of vector of cell groups and vector of cell names are different.\n");
 
 ret.names() = colnames;
  
 return ret;
}

/***********************************************************
 * 
 * Functions related with Sce format import/export
 *
 **********************************************************/
 
// Underlying internal function 

template <typename T>
void SceDataToBinMat(std::string fname,std::string ctype,bool isfull,bool transpose,Rcpp::NumericMatrix &MS4,
                     Rcpp::StringVector rownames, Rcpp::StringVector colnames,std::string comment)
{
 
 if (DEB & DEBSC)
 {
  Rcpp::Rcout << "Filling the internal " << (isfull ? "full" : "sparse") << " matrix...\n";
  Rcpp::Rcout.flush();
 }
 
 if (isfull)
 {
  FullMatrix<T> M(MS4.ncol(),MS4.nrow());
  
  if (transpose)
  {    
   for (indextype i=0; i<M.GetNRows(); i++)
    for (indextype j=0; j<M.GetNCols(); j++)
     M.Set(i,j,T(MS4(j,i)));
     
   if (ctype!="raw")
    M.SelfRowNorm(ctype);      // Normalization is done by rows since this matrix is really transposed from it origin
    
   if (DEB & DEBSC)
    Rcpp::Rcout << "Attaching vector of " << colnames.size() << " as row names and vector of " << rownames.size() << " as column names.\n";
    
   M.SetColNames(rownames);   
   M.SetRowNames(colnames);
  }
  else
  {
   for (indextype i=0; i<M.GetNRows(); i++)
    for (indextype j=0; j<M.GetNCols(); j++)
     M.Set(i,j,T(MS4(i,j)));
  
   if (ctype!="raw")
    M.SelfColNorm(ctype);      // Normalization is done as this functions says it does: by columns.
   
   if (DEB & DEBSC)
    Rcpp::Rcout << "Attaching vector of " << rownames.size() << " as row names and vector of " << colnames.size() << " as column names.\n";
   
   M.SetRowNames(rownames);  
   M.SetColNames(colnames);  
  } 
  
  if (comment!="")
   M.SetComment(comment);
 
  M.WriteBin(fname);
 } 
 else
 {
  SparseMatrix<T> M(MS4.ncol(),MS4.nrow());
  
  if (transpose)
  {    
   for (indextype i=0; i<M.GetNRows(); i++)
    for (indextype j=0; j<M.GetNCols(); j++)
     M.Set(i,j,T(MS4(j,i)));
     
   if (ctype!="raw")
    M.SelfRowNorm(ctype);      // Normalization is done by rows since this matrix is really transposed from it origin
    
   if (DEB & DEBSC)
    Rcpp::Rcout << "Attaching vector of " << rownames.size() << " as column names and vector of " << colnames.size() << " as row names.\n";
    
   M.SetColNames(rownames);   
   M.SetRowNames(colnames);
  }
  else
  {
   for (indextype i=0; i<M.GetNRows(); i++)
    for (indextype j=0; j<M.GetNCols(); j++)
     M.Set(i,j,T(MS4(i,j)));
  
   if (ctype!="raw")
    M.SelfColNorm(ctype);      // Normalization is done as this functions says it does: by columns.
   
   if (DEB & DEBSC)
    Rcpp::Rcout << "Attaching vector of " << colnames.size() << " as column names and vector of " << rownames.size() << " as row names.\n";
   
   M.SetRowNames(rownames);  
   M.SetColNames(colnames);  
  } 
  
  if (comment!="")
   M.SetComment(comment);
 
  M.WriteBin(fname);
 } 
}

//' SceToJMat
//'
//' Gets a numeric matrix of counts in the single cell experiment (sce) format and writes it to a disk file in the jmatrix binary format.\cr
//' To use this function you will have to extract yourself the matrix of counts (and may be the vectors of row names and column names) from the sce
//' or other object type. Plase, see the Details section
//'
//' The package BiocGenerics offers a facility to get the counts matrix, the function `counts, so usually you may load this package and use
//' counts(your_sce_object) as first argument. But sometimes not, and for example in the DuoClustering, you have\cr 
//' \cr
//' M<-your_object@assays$data@listData$counts\cr
//' \cr
//' to extract the counts matrix but in splatter you would have\cr
//' \cr
//' M<-your_object@assays@data@listData$counts\cr
//' \cr
//' (which is not exactly the same...)\cr
//'
//' The message, unfortunately, is: extract the data inspecting the internal structure of the object in the package that provided the data you are using.\cr
//' We assume, nevertheless, that if the matrix is M,\cr
//' attr(M,"dim")[1] is the number of rows (genes)\cr
//' attr(M,"dim")[2] is the number of columns (cells)\cr
//' attr(M,"dimnames")[1] is the vector of row names (names of genes)\cr
//' attr(M,"dimnames")[2] is the vector of colums names (names of cells)\cr
//' But if the matrix has not row or column names, or even if it has but you want to overwrite them, you can pass a value for parameter rownames or colnames
//' that will be honored. If you do not pass one or both the function will try to get them from the matrix attributes, as stated before. If they do not exist
//' as attributes in the matrix, they will be left empty.
//' 
//' The parameter transpose has the default value of FALSE. But don't forget to set it to TRUE if you want the cells
//' (which in single cell common practice are by columns) to be written by rows. This will be needed later to calculate
//' the dissimilarity matrix, if this is the next step of your workflow. See help of CalcAndWriteDissimilarityMatrix.
//'
//' @param M         The numeric matrix (extracted from the sce object as counts(theobject) or otherwise directly from the sce object).
//' @param fname     A string with the name of the binary output file
//' @param rownames  The vector of strings with the row names (extracted from the sce object, or set by the user). Default: empty vector (column names will be extracted from the matrix dimnames, if present)
//' @param colnames  The vector of strings with the column names (extracted from the sce object, or set by the user). Default: empty vector (row names will be extracted from the matrix dimnames, if present)
//' @param mtype     A string to indicate the matrix type: 'full' or 'sparse'. Default: 'sparse'
//' @param ctype     The string 'raw' or 'log1' to write raw counts or log(counts+1), or the normalized versions, 'rawn' and 'log1n', which normalize ALWAYS BY COLUMNS (before transposition, if requested to transpose). Default: raw
//' @param valuetype The data type to store the matrix. It must be one of the strings 'uint32', 'float' or 'double'. Default: float
//' @param transpose Boolean to indicate if the matrix should be transposed before writing. See Details for a comment about this. Default: FALSE
//' @param comment   A comment to be stored with the matrix. Default: "" (no comment)
//' @return   No return value, called for side effects (creates a file)
//' @examples
//' # Sorry, we cannot provide an example here, since it would need the load of the splatter package.
//' # Please, see the vignette for examples
//' @export
// [[Rcpp::export]]
void SceToJMat(Rcpp::NumericMatrix &M, std::string fname, 
		Rcpp::Nullable<Rcpp::StringVector> rownames=R_NilValue, Rcpp::Nullable<Rcpp::StringVector> colnames=R_NilValue,
		std::string mtype = "sparse", std::string ctype = "raw", std::string valuetype = "float", bool transpose = false, std::string comment="")
{
 if ((ctype != "raw") && (ctype != "log1") && (ctype != "rawn") && (ctype != "log1n"))
 {
  Rcpp::stop("The ctype argument can take only one of the string values 'raw', 'log1', 'rawn' and 'log1n'\n");
  return;
 }
 if ((mtype != "full") && (mtype != "sparse"))
 {
  Rcpp::stop("The mtype argument can take only one of the string values 'full' or 'sparse'\n");
  return;
 }
 bool isfull = (mtype == "full");
 
 if ((valuetype != "float") && (valuetype != "double") && (valuetype != "uint32"))
 {
  Rcpp::stop("The valuetype argument can take only one of the string values 'uint32', 'float' or 'double'\n");
  return;
 }
 if ((valuetype == "uint32") && ((ctype == "log1") || (ctype == "log1n")) )
 {
  Rcpp::stop("Rescaling as log(counts+1) requires output type to be float or double, not uint32, since it produces decimal numbers.\n");
  return;
 }
 if ((valuetype == "uint32") && (ctype == "rawn") )
 {
  Rcpp::stop("Rescaling as normalized raw counts requires output type to be float or double, not uint32, since it produces decimal numbers.\n");
  return;
 }
 
 unsigned int nrows=M.nrow();
 unsigned int ncols=M.ncol();
 if (DEB & DEBSC)
  Rcpp::Rcout << "The matrix of counts is of dimension [" << nrows << " x " << ncols << "].\n";
  
 if (nrows<2 || ncols<2)
 {
  Rcpp::stop("Count matrix must have at least two rows and two columns.\n");
  return;
 }
  
 Rcpp::StringVector rnames,cnames;
 
 if (M.hasAttribute("dimnames"))
 {
  Rcpp::List q4 = M.attr("dimnames");
  
  if (rownames.isNull())
  {
   rnames=q4[0];
   if (nrows!=(unsigned int)rnames.length())
    Rcpp::stop("Strange Matrix object. The number of rows in the matrix differs from the length of the vector of row names.\n");
   if (DEB & DEBSC)
    Rcpp::Rcout << "The passed matrix has row names and they will be used.\n";
  }
  else
  {
   if (DEB & DEBSC)
    Rcpp::Rcout << "The passed matrix had row names, but they will not be taken into accounts, since you have given a value to parameter rownames,\n";
   rnames=rownames;
  }
  
  if (colnames.isNull())
  {
   cnames=q4[1];
   if (ncols!=(unsigned int)cnames.length())
    Rcpp::stop("Strange Matrix object. The number of columns in the matrix differs from the length of the vector of column names.\n"); 
   if (DEB & DEBSC)
    Rcpp::Rcout << "The passed matrix has column names and they will be used.\n";
  }
  else
  {
   if (DEB & DEBSC)
    Rcpp::Rcout << "The passed matrix had column names, but they will not be taken into accounts, since you have given a value to parameter colnames,\n";
   cnames=colnames;
  }
 }
 else  // It the passed matrix has no row and column names, they will be filled with the passed values, if any.
 {
  rnames = (rownames.isNull()) ? Rcpp::StringVector() : rownames; 
  cnames = (colnames.isNull()) ? Rcpp::StringVector() : colnames;
 }
 
 if (valuetype=="uint32")
  SceDataToBinMat<unsigned int>(fname,ctype,isfull,transpose,M,rnames,cnames,comment);
 
 if (valuetype=="float")
  SceDataToBinMat<float>(fname,ctype,isfull,transpose,M,rnames,cnames,comment);
  
 if (valuetype=="double") 
  SceDataToBinMat<double>(fname,ctype,isfull,transpose,M,rnames,cnames,comment); 
}

/*
//' LogAndNormalize
//'
//' Takes the names of a file containing a binary matrix in jmatrix format and a file that will contain the transformed/normalized matrix, also in binary jmatrix format
//' The only possible transformation is the substitutiong of each value c by log2(c+1). AFTER that, normalization can be applied, either by rows or by columns.
//' Input matrix must be of type 'uint32' (unsigned long), 'float' or 'double'.
//' The transformed/normalized matrix will have the same character (full or sparse) of the original matrix and the datatype requested by the user.
//' Setting log1n to FALSE, nomralization to 'none' and otype to a type different from that of the input matrix can be used to change the data type.
//'
//' A call to this function with log1n=FALSE, normtype='none' and the ouput type equal to that of the input matrix would not alter the original matrix, and therefore
//' is not executed, raising an error.
//' Degrading the input matrix to a lower type (i.e.: 'double' to 'float' or to 'uint32' or 'float' to 'uint32' is allowed, but will provoke a precision loss and it raises a warning.
//' Apply the log1n transformation requires the output type be 'float' or 'double'. Otherwise, an error is raised.
//' This function cannot be applied to symmetric matrices. If ifname contains such a matrix, it would raise an error.
//'
//' @param ifname Name of the file containing the original matrix
//' @param ofname Name of the file containing the transformed matrix
//' @param log1n  Boolean. TRUE to apply the log2(c+1) transformation to any c value, FALSE to write the value c as is
//' @param normtype String to ask for the normalization to be applied with value 'none', 'byrows' or 'bycols'
//' @param otype  Data type of the output matrix. It can be 'uint32', 'float' or 'double'
//' @return   No return value, called for side effects (creates a file)
//' @examples
//' # Sorry, we cannot provide an example here, since it would need the load of the splatter package.
//' # Please, see the vignette for examples
//' @export
// [[Rcpp::export]]
void LogAndNormalize(std::string ifile, std::string ofile, bool log1n=false, std::string normtype, std::string otype)
{
 if ((normtype != "none") && (normtype != "byrows") && (normtype != "bycols"))
 {
  Rcpp::stop("normtype parameter must be 'none', 'byrows' or 'bycols'.\n");
  return;
 }
 if ((otype != "uint32") && (otype != "float") && (otype != "double"))
 {
  Rcpp::stop("otype parameter must be 'uint32', 'int32', 'float' or 'double'.\n");
  return;
 }
 unsigned char octype;
 if (otype=="uint32")
  octype=ULTYPE;
 if (otype=="float")
  octype=FTYPE;
 if (otype=="double")
  octype=DTYPE;
  
 unsigned char mtype,ctype,endian,mdinf;
 indextype nrows,ncols;
 
 MatrixType(ifile,mtype,ctype,endian,mdinf,nrows,ncols);
 
 if (mtype==MTYPESYMMETRIC)
 {
  Rcpp::stop("Input matrix cannot be symmetric, only full or sparse.\n");
  return;
 }
 
 if ((ctype!=ULTYPE) && (ctype!=FTYPE) && (ctype!=DTYPE))
 {
  Rcpp::stop("Input matrix is not of one of the allowed types (uint32, float or double).\n");
  return;
 }
 
 if ((log1n==true) && (otype=="uint32"))
 {
  Rcpp::stop("Cannot apply log1n transformation and write the result as integer. It should be float or double.\n");
  return;
 }
 
 if (
     (otype=="uint32" && (ctype=="float" || (ctype=="double"))) 
     ||
     (otype=="float") && (ctype=="double")
    )
    Rcpp::warn("Warning: you are downgrading data type. This will provoke a loss of precision.\n");
    
 if ((normtype=="none") && (log1n==false) && (ctype==octype))
 {
  Rcpp::stop("Not normalizing, not applying the log transformation and keeping the same datatype would not change the original file so it is not executed. In this case, just make a copy the file.\n");
  return;
 }
 
 if (mtype==MTYPEFULL)
 {
     switch (ctype)
     {
        case ULTYPE: { FullMatrix<unsigned long> M(ifile); WriteFullNorm(M,normtype,octype,ofile); break; };
        case FTYPE:  { FullMatrix<float> M(ifile);  WriteFullNorm(M,normtype,octype,ofile); break; };
        case DTYPE:  { FullMatrix<double> M(ifile); WriteFullNorm(M,normtype,octype,ofile); break; };
        default: break;
    }
 }
 if (mtype==MTYPESPARSE)
 {
     switch (ctype)
     {
        case ULTYPE: { SparseMatrix<unsigned long> M(ifile); WriteSparseNorm(M,normtype,octype,ofile); break; };
        case FTYPE:  { SparseMatrix<float> M(ifile);  WriteSparseNorm(M,normtype,octype,ofile); break; }; 
        case DTYPE:  { SparseMatrix<double> M(ifile); WriteSparseNorm(M,normtype,octype,ofile); break; };
        default: break;
    }
 }
 
}
*/
