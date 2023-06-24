#include <matmetadata.h>

extern unsigned char DEB;

/***********************************************************
 * 
 * Functions related with jmatrix metadata manipulation
 *
 **********************************************************/
 
// Internal functions to get metadata: names of rows and columns

// Auxiliary functions to skip the mark of binary block separation
unsigned char ChSep(std::ifstream &ifile)
{
 unsigned char dummy[BLOCKSEP_LEN];
 
 ifile.read((char *)dummy,BLOCKSEP_LEN);
 
 size_t i=0;
 while (i<BLOCKSEP_LEN && dummy[i]==BLOCKSEP[i])
  i++;
 
 return ((i==BLOCKSEP_LEN) ? READ_OK : ERROR_READING_SEP_MARK);
}

// Auxiliary functions to read metadata: row or column names (one after the other)
unsigned char RNames(std::ifstream &ifile,std::vector<std::string> &names)
{
 char dummy[MAX_LEN_NAME];
 char b;
 
 indextype j=0;
 do
 {
  b=char(ifile.get());
  
  if (ifile.eof())    // We have reached the end-of-file
   return ((j==0) ? READ_OK : ERROR_READING_STRINGS);
   
  if ((unsigned char)b==BLOCK_MARK)     // We have at the start between block of metadata
  {
   ifile.unget();
   return READ_OK;
  }
  
  if (b==0x00)   // We have reached end of this name
  {
   dummy[j]=b;   // Set it as end of string
   names.push_back(std::string(dummy));     // and take note of the string
   j=0;                             // Reinitialize the position in dummy
  }
  else
  { 
   dummy[j]=b;         // Usual case: still reading characters
   j++;
   if (j>=MAX_LEN_NAME)
    return ERROR_READING_STRINGS;
  }
 } 
 while (!ifile.eof());
 
 return ERROR_READING_STRINGS;
}

// Internal function to get simultaneously the row and column names
void InternalGetBinNames(std::string fname,unsigned char whichnames,std::vector<std::string> &rnames,std::vector<std::string> &cnames)
{
 unsigned char mtype,ctype,endian,mdinfo;
 indextype nrows,ncols;
 
 MatrixType(fname,mtype,ctype,endian,mdinfo,nrows,ncols);
 
 if (
     ( ((mdinfo & ROW_NAMES)==0) && ((whichnames & ROW_NAMES)!=0) ) ||
     ( ((mdinfo & COL_NAMES)==0) && ((whichnames & COL_NAMES)!=0) )
    )
 {
  if (DEB & DEBJM)
  {
   std::string w;
   if ((whichnames & ROW_NAMES) && (whichnames & COL_NAMES))
    w="Asking for row and colum names in file "+fname+", which did not store at least one of such data (even if there is one, the returned value will be empty).\n";
   else
   {
    if (whichnames & ROW_NAMES)
     w="Asking for row names in file "+fname+", which did not store such data.\n";
    else
     w="Asking for column names in file "+fname+", which did not store such data.\n";
   } 
   Rcpp::warning(w);
  }
  return;
 }
  
 unsigned long long start_metadata,start_comment;
 PositionsInFile(fname,&start_metadata,&start_comment);
 
 // Then, position the stream at the beginning of the metadata
 std::ifstream f(fname.c_str());
 f.seekg(start_metadata,std::ios::beg);
 if (whichnames & ROW_NAMES)
 {
  if (RNames(f,rnames) == ERROR_READING_STRINGS)
  {
   f.close();
   Rcpp::stop("Cannot read row names from binary file (even they are supposed to be there...).\n");
  }
  if (ChSep(f)==ERROR_READING_SEP_MARK)
   Rcpp::stop("Cannot read separation mark from binary file (even it should be supposed to be there...).\n");
 }
 else
 {
  // Row names are before column names in binary file. This obligues us to read them, if they are there, even if the caller did not ask for them.
  // But in this case parameter rnames will not be altered
  if (mdinfo & ROW_NAMES)
  {
   std::vector<std::string> dummy;
   if (RNames(f,dummy) == ERROR_READING_STRINGS)
   {
    f.close();
    Rcpp::stop("Cannot read row names from binary file (even they are supposed to be there...).\n");
   }
   if (ChSep(f)==ERROR_READING_SEP_MARK)
    Rcpp::stop("Cannot read separation mark from binary file (even it should be supposed to be there...).\n");
  }
 }
 
 if (whichnames & COL_NAMES)
 {
  if (RNames(f,cnames) == ERROR_READING_STRINGS)
  {
   f.close();
   Rcpp::stop("Cannot read column names from binary file (even they are supposed to be there...).\n");
  }
 }
 // On the contrary, if no column names have been asked, we don't need to read them.

 f.close();
}

//' GetJRowNames
//'
//' Returns a R StringVector with the row names of a matrix stored in the binary format of package jmatrix, if it has them stored.
//'
//' @param  fname  String with the file name that contains the binary data.
//' @return A R StringVector with the row names, or the empty vector if the binary file has no row names as metadata.
//' @examples
//' Rf <- matrix(runif(48),nrow=6)
//' rownames(Rf) <- c("A","B","C","D","E","F")
//' colnames(Rf) <- c("a","b","c","d","e","f","g","h")
//' tmpfile1=paste0(tempdir(),"/Rfullfloat.bin")
//' JWriteBin(Rf,tmpfile1,dtype="float",dmtype="full",comment="Full matrix of floats")
//' rn<-GetJRowNames(tmpfile1)
//' rn
//' @export
// [[Rcpp::export]] 
Rcpp::StringVector GetJRowNames(std::string fname)
{
 std::vector<std::string> rn,cn;
 InternalGetBinNames(fname,ROW_NAMES,rn,cn);
 Rcpp::StringVector rnames(rn.size());
 for (size_t i=0;i<rn.size();i++)
  rnames[i]=Rcpp::String(rn[i]);
 return rnames;
}

//' GetJColNames
//'
//' Returns a R StringVector with the column names of a matrix stored in the binary format of package jmatrix, if it has them stored.
//'
//' @param fname  String with the file name that contains the binary data.
//' @return A R StringVector with the column names, or the empty vector if the binaryfile has no column names as metadata.
//' @examples
//' Rf <- matrix(runif(48),nrow=6)
//' rownames(Rf) <- c("A","B","C","D","E","F")
//' colnames(Rf) <- c("a","b","c","d","e","f","g","h")
//' tmpfile1=paste0(tempdir(),"/Rfullfloat.bin")
//' JWriteBin(Rf,tmpfile1,dtype="float",dmtype="full",comment="Full matrix of floats")
//' cn<-GetJColNames(tmpfile1)
//' cn
//' @export
// [[Rcpp::export]] 
Rcpp::StringVector GetJColNames(std::string fname)
{
 std::vector<std::string> rn,cn;
 InternalGetBinNames(fname,COL_NAMES,rn,cn);
 Rcpp::StringVector cnames(cn.size());
 for (size_t i=0;i<cn.size();i++)
  cnames[i]=Rcpp::String(cn[i]);
 return cnames;
}

//' GetJNames
//'
//' Returns a R list of two elements, rownames and colnames, each of them being a R StringVector with the corresponding names
//'
//' @param fname  String with the file name that contains the binary data.
//' @return N["rownames","colnames"]: A list with two elements named rownames and colnames which are R StringVectors.
//'         If the binary file has no row or column names as metadata BOTH will be returned as empty vectors, even if one of them exists.
//'         If you want to extract only one, use either GetJRowNames or GetJColNames, as appropriate.
//' @examples
//' Rf <- matrix(runif(48),nrow=6)
//' rownames(Rf) <- c("A","B","C","D","E","F")
//' colnames(Rf) <- c("a","b","c","d","e","f","g","h")
//' tmpfile1=paste0(tempdir(),"/Rfullfloat.bin")
//' JWriteBin(Rf,tmpfile1,dtype="float",dmtype="full",comment="Full matrix of floats")
//' N<-GetJNames(tmpfile1)
//' N["rownames"]
//' N["colnames"]
//' @export
// [[Rcpp::export]]
Rcpp::List GetJNames(std::string fname)
{
 std::vector<std::string> rn,cn;
 InternalGetBinNames(fname,ROW_NAMES | COL_NAMES,rn,cn);
 
 Rcpp::StringVector rnames(rn.size());
 for (size_t i=0;i<rn.size();i++)
  rnames[i]=Rcpp::String(rn[i]);
  
 Rcpp::StringVector cnames(cn.size());
 for (size_t i=0;i<cn.size();i++)
  cnames[i]=Rcpp::String(cn[i]);

 Rcpp::List ret;
 ret["rownames"] = rnames;
 ret["colnames"] = cnames;

 return ret;
}

