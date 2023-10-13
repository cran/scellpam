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
 * The set of functions to import/export data form/to the most used formats in single cell:
 * Seurat/dgCMatrix and Sce (Single cell experiment).
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

