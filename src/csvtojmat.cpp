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
 * The set of functions to import/export data form/to CSV (comma separated values).
 * CSV does not need anything else, since it is self-described.
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
void CsvDataToBinMat(std::string ifname,std::string ofname,unsigned char vtype,std::string ctype,char csep,unsigned char vmtype,bool transpose,std::string comment)
{
 switch (vmtype)
 {
  case MTYPEFULL:
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
   break;
  }
  case MTYPESPARSE:
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
   break;
  }
  case MTYPESYMMETRIC:
  {
   SymmetricMatrix<T> M(ifname,vtype,csep);
   if (comment!="")
    M.SetComment(comment);
   M.WriteBin(ofname);
   break;
  }
  default: break;
 }
}

//' CsvToJMat
//'
//' Gets a csv/tsv file and writes to a disk file the binary matrix of counts contained in it in the jmatrix binary format.\cr
//' First line of the .csv is supposed to have the field names.\cr
//' First column of each line is supposed to have the row name.\cr
//' The fields are supposed to be separated by one occurrence of a character-field sepparator (usually, comma or tab)
//' .tsv files can be read with this function, too, setting the csep argument to '\\t'
//'
//' The parameter transpose has the default value of FALSE. But don't forget to set it to TRUE if you want the cells
//' (which in single cell common practice are by columns) to be written by rows. This will be needed later to calculate
//' the dissimilarity matrix, if this is the next step of your workflow. See help of CalcAndWriteDissimilarityMatrix
//' 
//' Special note for loading symmetric matrices:\cr
//' If you use this function to load what you expect to be a symmetric matrix from a .csv file, remember that the input table
//' MUST be square, but only the lower-diagonal matrix will be stored, including the main diagonal. The rest of the input table is
//' completely ignored, except to check that there are values in it. It is not checked if the table really represents a
//' symmetric matrix or not.\cr
//' Furthermore, symmetric matrices can only be loaded in raw mode, i.e.: no normalization is allowed, and they cannot be transposed.
//' 
//' @param ifname    A string with the name of the .csv/.tsv text file.
//' @param ofname    A string with the name of the binary output file.
//' @param mtype     A string to indicate the matrix type: 'full', 'sparse' or 'symmetric'. Default: 'sparse'
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
 if ((mtype != "full") && (mtype != "sparse") && (mtype != "symmetric"))
 {
  Rcpp::stop("The mtype argument can take only one of the string values 'full' or 'sparse'\n");
  return;
 }
 unsigned char vmtype = (mtype=="full") ? MTYPEFULL : ((mtype=="sparse") ? MTYPESPARSE : MTYPESYMMETRIC);
 if ((vmtype==MTYPESYMMETRIC) && (ctype!="raw"))
 {
  Rcpp::stop("Symmetric matrices cannot be normalized. Its ctype parameter must be left by default as 'raw'.");
  return;
 }
 if ((vmtype==MTYPESYMMETRIC) && (transpose==true))
 {
  Rcpp::stop("Symmetric matrices cannot be asked to be transposed (which is absurd, anyway, since they would remain invariant).\nBe sure your actual data are in the lower-diagonal part plus the main diagonal of the .csv table.");
  return;
 }
 
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
  CsvDataToBinMat<unsigned int>(ifname,ofname,ULTYPE,ctype,csep,vmtype,transpose,comment);
  
 if (valuetype=="float")
  CsvDataToBinMat<float>(ifname,ofname,FTYPE,ctype,csep,vmtype,transpose,comment);
  
 if (valuetype=="double")
  CsvDataToBinMat<double>(ifname,ofname,DTYPE,ctype,csep,vmtype,transpose,comment);
}

//' JMatToCsv
//'
//' Writes a binary matrix in the jmatrix package format as a .csv file. This is mainly for checking/inspection and
//' to load the data from R as read.csv, if the memory of having all data as doubles allows doing such thing.
//' 
//' The numbers are written to text with as many decimal places as allowed by its data type (internally obtained
//' with std::numeric_limits<type>::max_digits10)\cr
//' NOTE ON READING FROM R: to read the .csv files exported by this function you MUST use the R function read.csv
//' (not read.table) AND set its argument row.names to 1, since we always write a first column with the row names,
//' even if the binary matrix does not store them; in this case they are simply "1","2",...
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

