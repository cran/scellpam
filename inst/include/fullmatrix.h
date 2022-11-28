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

#ifndef FULLMATRIX_H
#define FULLMATRIX_H

#include <jmatrix.h>

/**
 * @FullMatrix class to hold full matrices (all space booked in memory)
 */
template <typename T>
class FullMatrix: public JMatrix<T>
{
 public:
    /**
     * Default constructor
     */
    FullMatrix();
    
    /**
     * Constructor with number of rows and columns
     * 
     * @param unsigned long nrows: number of rows
     * @param unsigned long ncols: number of columns
     */
    FullMatrix(indextype nrows,indextype ncols);

    /**
     * Copy constructor
     *
     * @param FullMatrix& other: Reference to the Matrix to be copied
     */
    FullMatrix( const FullMatrix<T>& other );

    /**
     * Constructor to fill the matrix content from a binary file
     * Binary file header as explained in the documetation to JMatrix::WriteBin
     * 
     * @param string fname: The name of the file to read
     * 
     */
    FullMatrix(std::string fname);
    
    /**
     * Constructor to fill the matrix content from a csv file
     * First line is supposed to have the field names and it is ignored.
     * First column of each line is supposed to have the field name, and it is ignored, too.
     * The passed character is the expected field sepparator (usually, comma or tab)
     * 
     * @param string fname: The name of the csv file to read
     * @param unsigned char vtype: The data type to be stored
     * @param char csep: The character used as field sepparator
     * 
     */
    FullMatrix(std::string fname,unsigned char vtype,char csep);
    
    /** 
     * Function to resize the matrix
     *  WARNING: previous content, if any, IS LOST (to be reviewed)
     * 
     * @param indextype newnr: new number of rows
     * @param indextype newnc: new number of cols
     * 
     */
    void Resize(indextype newnr,indextype newnc);
    
    /**
     * Destructor
     */
    ~FullMatrix();

    /**
     * Assignment operator
     *
     * @param FullMatrix& other: Reference to the Matrix to be assigned
     * @return Reference to the newly created Matrix
     */
    FullMatrix<T>& operator= ( const FullMatrix<T>& other );

    /**
     * Transpose-assignment
     * 
     * @param FullMatrix& other: Reference to the Matrix to be assigned
     * @return Reference to the newly created Matrix, which is the transpose of the passed one
     */
    FullMatrix<T>& operator!= ( const FullMatrix<T>& other);
    
    /** 
     * Function to get acess to an element
     * 
     * @param indextype r: The row to access
     * @param indextype c: The columns to access
     * 
     * @return T: value at (r,c) of matrix
     */
    T Get(indextype r,indextype c);
    
    /** 
     * Function to set an element
     * 
     * @param indextype r: The row to access
     * @param indextype c: The columns to access
     * @param T: the value to be set
     * 
     */
    void Set(indextype r,indextype c,T v);
    
    /**
      * Function to get a row as a pointer to the content type.
      * The pointer is not returned since that way it does not need to be booked. The pointer to hold result is passed as parameter
      * and it is supposed to be properly allocated.
      *
      * @param indextype r: The row to get
      * @param T *v: pointer to the result
      * 
      */
     void GetRow(indextype r,T *v);
     
    /**
      * Function to get a full row as a pointer to the content type plus a pointer to an array of marks.
      * The content will have the values. The array of marks will be changed OR'ing the passed mark to each place where there is a value.
      * The pointers to the values and marks are not returned since that way they do not need to be booked. They are passed as parameters
      * and both are supposed to be properly allocated.
      * This strange way of storing data has been choosen because it will be specially suitable to calculate distance between vectors
      *
      * @param indextype r: The row to get
      * @param T *v: pointer to the values
      * @param unsigned char *m: pointer to the marks
      * @param unsigned char s: the value to be OR'ed to each place at the mark array
      * 
      */
     void GetFullRow(indextype r,unsigned char *m,unsigned char s,T *v);
     
     /**
      * Function to get from a full row a pointer to an array of marks signalling where the non-zero elements are.
      * The array of marks will be changed OR'ing the passed mark to each place where there is a value.
      * The pointer to the marks is not returned since that way it do not need to be booked. It is passed as a parameter
      * and is supposed to be properly allocated.
      *
      * @param indextype r: The row to get
      * @param unsigned char *m: pointer to the marks
      * @param unsigned char s: the value to be OR'ed to each place at the mark array
      * 
      */
     void GetMarksOfFullRow(indextype r,unsigned char *m,unsigned char s);
    
     /**
      * Function to alter the internal values of the matrix so that each row is normalized according to the requested normalization type
      * The purpose of this function can be achieved with a loop using Set and Get, but using the internal structure makes the task much faster
      *
      * @param string ctype:    The requested type of normalization: rawn, log1 or log1n
      *
      */
     void SelfRowNorm(std::string ctype);
     
     /**
      * Function to alter the internal values of the matrix so that each column is normalized according to the requested normalization type
      * The purpose of this function can be achieved with a loop using Set and Get, but using the internal structure makes the task much faster
      *
      * @param string ctype:    The requested type of normalization: rawn, log1 or log1n
      *
      */
     void SelfColNorm(std::string ctype);
       
    /**
     * Function to write the matrix content to a CSV file
     * 
     *  @param string fname:    The name of the file to write
     *  @param char csep:       The separator character between fields (default: , (comma))
     *  @param withquotes:      bool to indicate if field names in .csv must be written surrounded by quotes.
     * 
     */
    void WriteCsv(std::string fname,char csep=',',bool withquotes=false);
    
    
    /**
     * Function to write the matrix content to a binary file
     * See format at documentation of JMatrix::WriteBin
     * 
     *  @param string fname: The name of the file to write
     * 
     */
    void WriteBin(std::string fname);
    
    /**
     * Function to get memory in MB used by this full matrix
     * 
     */
    float GetUsedMemoryMB();
    
 private:
     T **data;
};

#endif // MATRIX_H
