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

#ifndef JMATRIX_H
#define JMATRIX_H

#include <sys/stat.h>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>		//std::remove_copy
#include <type_traits>
#include <Rcpp.h>
#include <debugpar.h>
#include <templatemacros.h>

///@{
/** 
 *       Constants for the possible matrix types
 *         Currently, they are no type (for errors), full matrix, sparse matrix and symmetric matrix.
 *
 */
const unsigned char MTYPENOTYPE=0x0F;		/*!< No matrix type */	
const unsigned char MTYPEFULL=0x00;		/*!< Full matrix */	
const unsigned char MTYPESPARSE=0x01;		/*!< Sparse matrix */	
const unsigned char MTYPESYMMETRIC=0x02;	/*!< Symmetric matrix */	
///@}

///@{
/** 
 *        Constants for the possible data types a matrix can hold.
 *         These are (apart of the no type for errors) integer types:
 *            char (8 bits), short int (16 bits), int (32 bits), long (32 bits), long long (64 bits) with their signed versions
 *         and float types in IEEE-754 format:
 *            float (32 bits), double (64 bits) and long double (128 bits)
 *
 */
const unsigned char NOTYPE=0x0F;	/*!< No data type */		
const unsigned char UCTYPE=0x00;	/*!< unsigned char */			
const unsigned char SCTYPE=0x01;	/*!< char */			
const unsigned char USTYPE=0x02;	/*!< unsigned short */			
const unsigned char SSTYPE=0x03;	/*!< short */			
const unsigned char UITYPE=0x04;	/*!< unsigned int */			
const unsigned char SITYPE=0x05;	/*!< int */			
const unsigned char ULTYPE=0x06;	/*!< unsigned long */			
const unsigned char SLTYPE=0x07;	/*!< long */			
const unsigned char ULLTYPE=0x08;	/*!< unsigned long long */			
const unsigned char SLLTYPE=0x09;	/*!< long long */			
const unsigned char FTYPE=0x0A;		/*!< float */			
const unsigned char DTYPE=0x0B;		/*!< double */				
const unsigned char LDTYPE=0x0C;	/*!< long double */	
///@}

///@{		
/**
 *       Constants for the possible endianness of a machine
 *         Only two values allowed for big and little endian
 *
 */
const unsigned char BIGEND=0x00;	/*!< Big endian*/			
const unsigned char LITEND=0xF0;	/*!< Little endian */			
///@}

///@{
/**
*	Constants for information about the metadata included with the matrix
*       Possible metadata stored are currently names of rows and names of columns
*       Errors are returned with these constants in case of bad reads
*
*/
const int READ_OK=0;
const int ERROR_READING_STRINGS=1;
const int ERROR_READING_ROW_NAMES=2;
const int ERROR_READING_COL_NAMES=3;
const int ERROR_READING_SEP_MARK=4;

const unsigned int  MAX_LEN_NAME=1023;
const unsigned int  COMMENT_SIZE=1024;      // All comments have fixed length. If no comment is used this will be a chunck of null characters.

const unsigned char NO_METADATA=0x00;
const unsigned char ROW_NAMES=0x01;
const unsigned char COL_NAMES=0x02;
const unsigned char COMMENT=0x04;

const unsigned int BLOCKSEP_LEN=4;
const unsigned char BLOCK_MARK=0xFF;
const unsigned char BLOCKSEP[BLOCKSEP_LEN]={BLOCK_MARK,0x45,0x42,BLOCK_MARK};  // The is 0xFF, E, B, 0xFF
///@}

/**
 * Returns the endianness of the machine where this function is called
 *
 * @return unsigned char: either the constant BIGEND or the constant LITEND
 */
unsigned char ThisMachineEndianness();

/**
 * Returns the file size of a file in a sufficiently large number (unsigned long long)
 * in a way (hopefully) independent of the operating system and of the architecture
 *
 * @param std::string File path
 * @return unsigned long long File size
 */
unsigned long long GetFileSize(std::string fname);

/**
 * Returns the positions of the start of metadata and start of comments (included inside metadata)
 * as absolute positions measured in bytes from the beginning of the file
 *
 * @param std::string File path
 * @param unsigned long long *start_of_metadata
 * @param unsigned long long *start_of_comment
 */
void PositionsInFile(std::string fname,unsigned long long *start_of_metadata,unsigned long long *start_of_comment);

/*! \brief Auxiliary functions to be used for error printing.
 *
 */

/**
 * Returns the name of the matrix whose type (as indentifier) is passed
 *
 * @param unsigned char: matrixtypeident. The type of the matrix (full, sparse, symmetric) as defined by the former constants.
 * @return string: A human-meaningful string decribing the type
 */
std::string MatrixTypeName(unsigned char matrixtypeident);

/**
 * Returns the name of the data type whose type (as indentifier) is passed
 *
 * @param unsigned char: datatypeident. The data type of the data in the matrix (unsigned char, char,...) as defined by the former constants.
 * @return string: A human-meaningful string decribing the data type
 */
std::string DataTypeName(unsigned char datatypeident);

/**
 * Returns a message to interpret the presence of metadata
 *
 * @param unsigned char: metadatainfo. The constant for information on which metadata are present as defined by the former constants.
 * @return string: A human-meaningful string decribing the present metadata
 */
std::string MetadataInfo(unsigned char metadatainfo);

/**
 * Returns the size in bytes of the data type whose type (as indentifier) is passed
 *
 * @param unsigned char: datatypeident. The data type of the data in the matrix (unsigned char, char,...) as defined by the former constants.
 * @return int: The size in bytes of one element of the passed data type.
 */
int SizeOfType(unsigned char datatypeident);

const unsigned short HEADER_SIZE=128;	/*!< The header size. We fix a header of 128 bytes. We don't need so much, but just in case in the future... */

/**
 * \typedef This is the type of the indexes. It is left fixed, but using a typedef allows easy recompilation of the library if
 * we decide otherwise in the future. Reducing to unsigned short would save memory, if you are sure that no matrix will never be
 * bigger than 65536, neither in rows nor in columns. Increasing it to 128 bits (unsigned long long) would be possible, but if you
 * need such size of matrices, better don't use this library...
 */
typedef unsigned int indextype;

/**
 * @JMatrix Wrapper class for all types of matrices. It is meant to hold some basic operations common to all of them
 * Even the instances of this class may now hold real data (the metadata, row and column names), they don't do until an "authentic" matrix is constructed.
 */
template <typename T>
class JMatrix
{
 public:
     /**
     * Default constructor
     *
     * @param unsigned_char_mtype: the matrix type (see constants at jmatrix.h)
     */
    JMatrix(unsigned char mtype);
    
    /**
     * Constructor with number of rows and columns
     * 
     * @param unsigned_char_mtype: the matrix type (see constants at jmatrix.h)
     * @param unsigned_long_nrows: number of rows
     * @param unsigned_long_ncols: number of columns
     */
    JMatrix(unsigned char mtype,indextype nrows,indextype ncols);

    /**
     * Constructor to fill the matrix content from a binary file
     * 
     * Binary file header as explained in the documentation to WriteBin
     * 
     * TODO PRELIMINARY VERSION. ASSUMES SAME ENDIANESS FOR WRITER AND READER MACHINE
     * 
     * @param string_fname: The name of the file to read
     * @param unsigned_char_mtype: the matrix type (see constants at jmatrix.h)
     * 
     */
    JMatrix(std::string fname,unsigned char mtype);
    
    /**
     * Constructor to fill the matrix content from a csv text file
     * 
     * 
     * @param string_fname: The name of the file to read
     * @param unsigned_char_mtype: the matrix type (see constants at jmatrix.h)
     * @param unsigned_char_vtype: the data type to be contained in the matrix
     * @param char_csep: The character expected to be the sepparator.
     * 
     */
    JMatrix(std::string fname,unsigned char mtype,unsigned char vtype,char csep);
    
    /**
     * Copy constructor
     *
     * @param JMatrix&_other: Reference to the Matrix to be copied
     */
    JMatrix(const JMatrix<T>& other );
    
    /**
     * Assignment operator
     *
     * @param JMatrix&_other: Reference to the Matrix to be assigned
     * @return JMatrix&: Reference to the newly created Matrix
     */
    JMatrix<T>& operator= ( const JMatrix<T>& other );

    /**
     * Transpose-assignment
     * 
     * @param JMatrix&_other: Reference to the Matrix to be assigned
     * @return JMatrix&: Reference to the newly created Matrix, which is the transpose of the passed one
     */
    JMatrix<T>& operator!= ( const JMatrix<T>& other);
    
    /**
     * Function to get number of rows
     * 
     * @return indextype: Number of rows of the matrix
     */
    indextype GetNRows() { return nr; };
    
    /**
     * Function to get number of columns
     * 
     * @return indextype: Number of columns of the matrix
     */
    indextype GetNCols() { return nc; };
    
    /**
     * Function to resize the matrix
     * WARNING: previous content, if any, IS LOST (TODO: to be reviewed)
     * 
     * @param indextype_newnr: new number of rows
     * @param indextype_newnc: new number of cols
     * 
     */
    void Resize(indextype newnr,indextype newnc);
    
    /** 
     * Function to get the matrix column names, if present
     *
     * @return: the vector of strings with the row names. Empty vector if not present
     *
     */
    std::vector<std::string> GetColNames();
    
     /** 
     * Function to get the matrix row names, if present
     *
     * @return: the vector of strings with the row names. Empty vector if not present
     *
     */
    std::vector<std::string> GetRowNames();
    
    /** 
     * Function to set the matrix column names.
     *
     * @param cnames: the std::vector of strings with the column names.
     *
     */
    void SetColNames(std::vector<std::string> cnames);
    
    /** 
     * Second function to set the matrix column names.
     *
     * @param cnames: the Rcpp::StringVector with the column names.
     *
     */
    void SetColNames(Rcpp::StringVector cnames);
    
     /** 
     * Function to set the matrix row names.
     *
     * @param rnames: the std::vector of strings with the row names.
     *
     */ 
    void SetRowNames(std::vector<std::string> rnames);
    
    /** 
     * Second function to set the matrix row names.
     *
     * @param rnames: the Rcpp::StringVector with the row names.
     *
     */
    void SetRowNames(Rcpp::StringVector rnames); 
    
    /** 
     * Function to get a string with the matrix comment, if any (or the empty string otherwise).
     *
     * @return: A string with the content of the comment area.
     *
     */
    std::string GetComment();
     
    /**
     * Function to set the matrix comment
     *
     * @param cm: A std::string with the comment to be stored.
     *
     */
    void SetComment(std::string cm);
    
    /**
     * Function to write the matrix content to a CSV  file
     * 
     *  @param string_fname: The name of the file to write
     *  @param char_csep:    The separator character between fields (default: , (comma))
     *  @param withquotes:   bool to indicate if field names in .csv must be written surrounded by quotes.
     * 
     */
    void WriteCsv(std::string fname, char csep=',',bool withquotes=false);
    
    /**
     * Function to write the matrix content to a binary file
     * 
     * The binary header will contain:
     *  - unsigned char t: matrix type (normal, sparse, symmetric). Other types might be added later
     *  - unsigned char dt: the data type of the matrix elements, which is one of
     *       unsigned/signed char, unsigned signed short, unsigned/signed long, unsigned/signed longlong, float,  double, longdouble 
     *    in its lower 4 bits and the endianess in its upper 4 bits (big or little).
     *  - indextype nr: number of rows
     *  - indextype nc: number of columns
     *  - unsigned char mdinfo: information on the presence/absence of metadata, currently row and/or column names.
     *
     *  This means that the size of the header is 2+2*sizeof(indextype)+1+empty_space.
     *  We have fixed the empty space so that total header size be 128 bytes.
     *
     *  Then comes the content as raw data, by rows.
     *  Then, the row names, as character arrays separated by null character (0x00)
     *  Then, a separation mark (the succesion of bytes 0xFF 0x45 0x42 0xFF). 0x45 and 0x42 are ASCII characters EB, for End Block
     *  Then, the column names, as character arrays separated by null character (0x00)
     * 
     *  HEADER AND DATA ARE ALWAYS WRITTEN IN THE ENDIANESS OF THE MACHINE WHICH EXECUTES THIS CODE
     *  Nevertheless, other machines will know about it, since it is declared in the first byte of the written file.
     * 
     *  @param string fname: The name of the file to write
     *  @param unsigned char mtype: The type identifier
     * 
     */
    void WriteBin(std::string fname,unsigned char mtype);
    
 protected: 
 	indextype nr,nc;
 	unsigned char jctype;
 	std::ifstream ifile;
 	std::ofstream ofile;
 	unsigned char TypeNameToId();
 	bool ProcessDataLineCsv(std::string line,char csep,T *rowofdata);
 	bool ProcessDataLineCsvForSymmetric(std::string line,char csep,indextype nrow,std::vector<T> &rowofdata);
 	int ReadMetadata();
 	void WriteMetadata();
 	std::vector<std::string> rownames;
 	std::vector<std::string> colnames;
 	char comment[COMMENT_SIZE];
 	void SetDataType(unsigned char dtype);
 	std::string CleanQuotes(std::string s);
 	
 private:
 	unsigned char jmtype;
 	unsigned char mdinfo;
 	bool ProcessFirstLineCsv(std::string line,char csep);
	void WriteNames(std::vector<std::string> &names);
	indextype ReadNames(std::vector<std::string> &names);
	indextype CheckSep();
};

/**
 *      Constant defined to check access to matrix elements
 *
 *      The simple fact of being defined at compilation time adds a test in each matrix access to be sure we are
 *      not out of bound; if we are, a run-time error is raised.
 *
 *      This is obviously safer but at the expense of adding overhead and a slight increment of run time.\n
 *      Comment this constant if you are absolutely sure your program does not make any Get or Set out of bounds.
 *
 */
//#define WITH_CHECKS_MATRIX
 
/**
 * Auxiliary funcions to check characteristics of matrix stored in binary file looking only at its header.
 *
 * @param string fname			Name of the binary file
 * @param unsigned char &mtype		Returns matrix type (full,sparse,symmetric)
 * @param unsigned char &ctype  	Returns matrix data type
 * @param unsigned char &endianness	Returns endianness (big,little) of the data stored in the matrix
 * @param indextype     &nrows		Returns the number of rows
 * @param indextype     &ncols		Returns the number of columns
 * @param metadatainfo  &mdinf          Returns the signal for the presence of metadata
 *
 */
void MatrixType(std::string fname,unsigned char &mtype); 
void MatrixType(std::string fname,unsigned char &mtype,unsigned char &ctype);
void MatrixType(std::string fname,unsigned char &mtype,unsigned char &ctype,unsigned char &endianness); 
void MatrixType(std::string fname,unsigned char &mtype,unsigned char &ctype,unsigned char &endianness,unsigned char &mdinf); 
void MatrixType(std::string fname,unsigned char &mtype,unsigned char &ctype,unsigned char &endianness,unsigned char &mdinf,indextype &nrows,indextype &ncols);

/** 
 * Auxiliary function to add or remove quotes as needeed
 * 
 * @param string s          String to fix if needed
 * @param bool   withquotes Boolean to indicate if the returned string must be surrounded by quotes or not
 * @return string           The string fixed as needed
 *
 */
 std::string FixQuotes(std::string s,bool withquotes);
#endif // JMATRIX_H
