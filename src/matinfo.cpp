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

/***********************************************************
 * 
 * Functions related with jmatrix binary format information
 *
 **********************************************************/
 
// ******** Function to get information of the matrix *************

//' JMatInfo
//'
//' Shows in the screen or writes to a file information about a matrix stored in the binary format of package jmatrix
//'
//' @param fname  String with the file name that contains the binary data.
//' @param fres   String with the name of the file to write the information. Default: "" (information is written to the console)
//' @return       No return value, called for its side effects (writes on screen or creates a file)
//' @examples
//' Rf <- matrix(runif(48),nrow=6)
//' rownames(Rf) <- c("A","B","C","D","E","F")
//' colnames(Rf) <- c("a","b","c","d","e","f","g","h")
//' tmpfile1=paste0(tempdir(),"/Rfullfloat.bin")
//' JWriteBin(Rf,tmpfile1,dtype="float",dmtype="full",comment="Full matrix of floats")
//' JMatInfo(tmpfile1)
//' @export
// [[Rcpp::export]]
void JMatInfo(std::string fname, std::string fres = "")
{
 unsigned char mtype,ctype,endian,mdinfo;
 indextype nrows,ncols;
 
 MatrixType(fname,mtype,ctype,endian,mdinfo,nrows,ncols);
 
 char comment[COMMENT_SIZE];
 unsigned long long start_metadata,start_comment;
 PositionsInFile(fname,&start_metadata,&start_comment);
 
 if (mdinfo & COMMENT)
 {
  std::ifstream f(fname.c_str());
  f.seekg(start_comment,std::ios::beg);
  f.read((char *)comment,COMMENT_SIZE);
  f.close();
 }
 
 std::streambuf *buf;
 std::ofstream fr;
 if (fres != "")
 {
  fr.open(fres.c_str());
  if (!fr.is_open())
  {
   std::ostringstream errst;
   errst << "File " << fres << " cannot be opened to write.\n";
   Rcpp::stop(errst.str());
   return;
  }
  buf = fr.rdbuf();
 }
 else
  buf = Rcpp::Rcout.rdbuf();
  
 std::ostream out(buf);
  
 out << "File:               " << fname << std::endl;
 out << "Matrix type:        ";
 switch (mtype)
    {
        case MTYPEFULL:         out << "FullMatrix\n"; break;
        case MTYPESPARSE:       out << "SparseMatrix\n"; break;
        case MTYPESYMMETRIC:    out << "SymmetricMatrix\n"; break;
        default:                out << "UnknownTypeMatrix\n"; break;
    }
 out << "Number of elements: " << (unsigned long)nrows*(unsigned long)ncols;
 if (mtype==MTYPESYMMETRIC)
  out << " (" << (unsigned long)nrows*(unsigned long)(ncols+1)/2 << " really stored)";
 out << std::endl;
 out << "Data type:          ";
 switch (ctype)
    {
        case UCTYPE: out <<  "unsigned char\n"; break;
        case SCTYPE: out <<  "char\n"; break;
        case USTYPE: out <<  "unsigned short int\n"; break;
        case SSTYPE: out <<  "short int\n"; break;
        case UITYPE: out <<  "unsigned int\n"; break;
        case SITYPE: out <<  "int\n"; break;
        case ULTYPE: out <<  "unsigned long\n"; break;
        case SLTYPE: out <<  "long\n"; break;
        case FTYPE:  out <<  "float\n"; break;
        case DTYPE:  out <<  "double\n"; break;
        case LDTYPE: out <<  "long double\n"; break;
        default: out <<  "unknown\n"; break;
    }
 out << "Endianness:         " << ((endian == BIGEND) ? "big endian" : "little endian");
 if (ThisMachineEndianness() == endian)
  out << " (same as this machine)\n";
 else
  out << " which is DIFFERENT from that of this machine.\n";
  
 out << "Number of rows:     " << nrows << std::endl;
 out << "Number of columns:  " << ncols << std::endl;
 out << "Metadata:           ";
 if (mdinfo==NO_METADATA)
  out << "None\n";
 else
 { 
  if ((mdinfo & ROW_NAMES) && (!(mdinfo & COL_NAMES)))
   out << "Stored only names of rows.\n";
  if (!(mdinfo & ROW_NAMES) && (mdinfo & COL_NAMES))
   out << "Stored only names of columns.\n";
  if ((mdinfo & ROW_NAMES) && (mdinfo & COL_NAMES))
   out << "Stored names of rows and columns.\n";
 }
 if (mdinfo & COMMENT)
  out << "Metadata comment:  \"" << comment << "\"\n";
  
 if (mtype==MTYPESPARSE)
 {
     unsigned long long full_size=(unsigned long long)nrows*(unsigned long long)ncols*SizeOfType(ctype);
     unsigned long long used_size=start_metadata-HEADER_SIZE;
     out << "Binary data size:   " <<  used_size << " bytes, which is " << 100.0*float(used_size)/float(full_size) << "% of the full matrix size (which would be " << full_size  << " bytes).\n";
 }
 
 if (fres != "")
  fr.close();
  
}

