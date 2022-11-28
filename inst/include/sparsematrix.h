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

#ifndef SPARSEMATRIX_H
#define SPARSEMATRIX_H

#include <jmatrix.h>

#include <algorithm>    // std::sort, std::stable_sort

template <typename T>
void sort_indexes_and_value(const std::vector<T> &v,std::vector<size_t> &idx,std::vector<indextype> &idv);

enum TrMark { transpose=0 };

/**
 * @SparseMatrix Class to hold arbitrarily big sparse matrices. Elements are stored with column index + value in a vector associated to each row.
 *               Time to set and get elements are of order O(log_2(Nc)) being Nc the number of columns.
 *               Space is O(N*(sizeof(element)+sizeof(index)), being element the type of the matrix contents and index that of the matrix index
 *               (which is currently unsigned int).
 */
template <typename T>
class SparseMatrix: public JMatrix<T>
{
 public:
    /**
     * Default constructor
     */
    SparseMatrix();

    /**
     * Constructor with number of rows and columns
     * 
     * @param indextype nrows: number of rows
     * @param indextype ncols: number of columns
     */
    SparseMatrix(indextype nrows,indextype ncols);
    
    /**
     * Constructor to fill the matrix content from a binary file
     * 
     * Binary file header as explained in the documetation to WriteBin
     * 
     * PRELIMINARY VERSION. ASSUMES SAME ENDIANESS FOR WRITER AND READER MACHINE
     * 
     * @param string fname: The name of the file to read
     * 
     */
    SparseMatrix<T>(std::string fname);
    
    /**
     * Constructor to fill the matrix content from the transpose of the matrix contained in a binary file
     * 
     * Binary file header as explained in the documentation to WriteBin
     * 
     * PRELIMINARY VERSION. ASSUMES SAME ENDIANESS FOR WRITER AND READER MACHINE
     * 
     * @param string fname: The name of the file to read
     * 
     */
    SparseMatrix<T>(std::string fname,TrMark);
    
    /**
     * Constructor to fill the matrix content from a csv file
     * First line is supposed to have the field names and is ignored.
     * First column of each line is supposed to have the field name, and is ignored, too.
     * The passed character is the expected field sepparator (usually, comma or tab)
     * 
     * @param string fname: The name of the csv file to read
      * @param unsigned char vtype: The data type to be stored
     * @param char csep: The character used as field sepparator
     * 
     */
    SparseMatrix(std::string fname,unsigned char vtype,char csep);
    
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
     * Copy constructor
     *
     * @param SparseMatrix& other: Reference to the SparseMatrix to be copied
     */
    SparseMatrix(const SparseMatrix& other);

    /**
     * Destructor
     */
    ~SparseMatrix();

    /**
     * Assignment operator
     *
     * @param SparseMatrix& other: Reference to the SparseMatrix to be assigned
     * @return Reference to the newly created SparseMatrix
     */
    SparseMatrix<T>& operator=(const SparseMatrix<T>& other);

    /**
     * Transpose-assignment
     * 
     * @param SparseMatrix& other: Reference to the SparseMatrix to be assigned
     * @return Reference to the newly created SparseMatrix, which is the transpose of the passed one
     */
    SparseMatrix<T>& operator!=(const SparseMatrix<T>& other);
    
    /** 
     * Function to get acess to an element
     * 
     * @param indextype r: The row to access
     * @param indextype c: The columns to access
     * 
     * @return T: value at (r,c) of matrix
     */
    T Get(indextype r,indextype c) const;
    
    /** 
     * Function to set an element
     * 
     * @param indextype r: The row of the element to be set
     * @param indextype c: The column of the element to be set
     * @param T: the value to be set
     * 
     */
    void Set(indextype r,indextype c,T v);
    
    /** 
     * Function to set a row (as two vectors of locations and values)
     * 
     * @param indextype r: The row to be set
     * @param vector<indextype> cv: The vector with the columns to be set
     * @param vector<T>: The vector with the corresponding values to be set 
     * 
     */
     void SetRow(indextype r,std::vector<indextype> vc,std::vector<T> v);
     
     /**
      * Function to get a row as a pointer to the content type. Row is not a sparse but a full vector with zeros when needed.
      * The pointer is not returned since that way it does not need to be booked. The pointer to hold result is passed as parameter
      * and it is supposed to be properly allocated.
      *
      * @param indextype r: The row to get
      * @param T *v: pointer to the result
      * 
      */
     void GetRow(indextype r,T *v);
    
     /**
      * Function to get a sparse row as a pointer to the content type plus a pointer to an array of marks.
      * The content will have the values. The array of marks will be changed OR'ing the passed mark to each place where there is a value.
      * The pointers to the values and marks are not returned since that way they do not need to be booked. They are passed as parameters
      * and both are supposed to be properly allocated.
      * This strange way of storing data has been choosen because it will be specially suitable to calculate distance between sparse vectors
      *
      * @param indextype r: The row to get
      * @param T *v: pointer to the values
      * @param unsigned char *m: pointer to the marks
      * @param unsigned char s: the value to be OR'ed to each place at the mark array
      * 
      */
     void GetSparseRow(indextype r,unsigned char *m,unsigned char s,T *v);
     
     /**
      * Function to get from a sparse row a pointer to an array of marks signalling where the non-zero elements are.
      * The array of marks will be changed OR'ing the passed mark to each place where there is a value.
      * The pointer to the marks is not returned since that way it do not need to be booked. It is passed as a parameter
      * and is supposed to be properly allocated.
      *
      * @param indextype r: The row to get
      * @param unsigned char *m: pointer to the marks
      * @param unsigned char s: the value to be OR'ed to each place at the mark array
      * 
      */
     void GetMarksOfSparseRow(indextype r,unsigned char *m,unsigned char s);
     
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
     * Function to write the matrix content to a CSV  file
     * 
     *  @param string fname: The name of the file to write
     *  @param char_csep:    The separator character between fields (default: , (comma))
     *  @param withquotes:      bool to indicate if field names in .csv must be written surrounded by quotes.
     * 
     */
     void WriteCsv(std::string fname,char csep=',',bool withquotes=false);
    
     /**
     * Function to write the matrix content into a binary file
     * For the header format, see the documentation of JMatrix
     *
     *  After the header comes the content as raw data, by rows, with this content for each row:
     *   - indextype ncr: number of non-zero entries of this row
     *   - ncr values of indextype with the numbers of the columns of this row occupied by non-zero entries
     *   - ncr elements of the current value type (the values of all non-zero entries of this row).
     * 
     *  @param string fname: The name of the file to write
     */
    void WriteBin(std::string fname);
     
    /**
     * Function to get memory in MB used by this sparse matrix (including values and additional indexes)
     * 
     */
    float GetUsedMemoryMB();
    
private:
    std::vector<std::vector<indextype>> datacols;
    std::vector<std::vector<T>> data;
};

#endif // SPARSEMATRIX_H
