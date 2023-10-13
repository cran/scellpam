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

#ifndef SYMMETRICMATRIX_H
#define SYMMETRICMATRIX_H

#include "jmatrix.h"

/**
 * @SymmetricMatrix Class to hold arbitrarily big symmetric square matrices. For a matrix of size NxN, only Nx(N+1)/2 elements are stored.
 */
template <typename T>
class SymmetricMatrix: public JMatrix<T>
{
 public:
    /**
     * Default constructor
     */
    SymmetricMatrix();

    /**
     * Constructor with number of rows/columns (the same, it is square)
     * 
     * @param indextype nrows: number of rows
     */
    SymmetricMatrix(indextype nrows);

    /**
     * TODO: make a copy constructor??
     */
     
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
    SymmetricMatrix(std::string fname); 
    
    /**
     * Constructor to fill the matrix content from a csv file
     * First line is supposed to have the field names and it is ignored.
     * First column of each line is supposed to have the field name, and it is ignored, too.
     * The passed character is the expected field sepparator (usually, comma or tab)
     * WARNING: even the .csv table must be square, only the lower-diagonal part (inclding main diagonal) is stored.
     * 
     * @param string fname: The name of the csv file to read
     * @param unsigned char vtype: The data type to be stored
     * @param char csep: The character used as field sepparator
     * 
     */
    SymmetricMatrix(std::string fname,unsigned char vtype,char csep);
    
    /** 
     * Function to resize the matrix
     *  WARNING: previous content, if any, IS LOST (to be reviewed)
     * 
     * @param indextype newnr: new number of rows (and columns)
     * 
     */
    void Resize(indextype newnr);
    
    /**
     * Copy constructor
     *
     * @param SymmetricMatrix& other: Reference to the SymmetricMatrix to be copied
     */
    SymmetricMatrix(const SymmetricMatrix<T>& other);

    /**
     * Destructor
     */
    ~SymmetricMatrix();

    /**
     * Assignment operator
     *
     * @param SymmetricMatrix& other: Reference to the SymmetricMatrix to be assigned
     * @return Reference to the newly created SymmetricMatrix
     */
    SymmetricMatrix<T>& operator=(const SymmetricMatrix<T>& other);

    /**
     * Test of correctness
     * This is meant to test if the symmetric matrix is a distance or dissimilarity matrix.
     * It checks that all elements in the main diagonal are 0 and all outside the main diagonal are stricty possitive
     * 
     * @return true if the matrix can be a distance or dissimilarity matrix. false otherwise.
     * 
     */
    bool TestDistDisMat();
    
    /** 
     * Function to get acess to an element
     * 
     * @param indextype r: The row to access
     * @param indextype c: The columns to access
     * 
     * @return T: value at (r,c) of matrix
     */
#ifdef WITH_CHECKS_MATRIX
    T Get(indextype r,indextype c);
#else
    inline T Get(indextype r,indextype c) { return (c<=r) ? data[r][c] : data[c][r]; };
#endif
    /** 
     * Function to set an element
     * 
     * @param indextype r: The row to access
     * @param indextype c: The columns to access
     * @param T: the value to be set
     * 
     */
#ifdef WITH_CHECKS_MATRIX
    void Set(indextype r,indextype c,T v); 
#else  
    inline void Set(indextype r,indextype c,T v) { if (c<=r) data[r][c]=v; else data[c][r]=v; };
#endif

    /**
     * Function to get the sum of a row (used frequently by PAM)
     * 
     * @param indextype r: The row whose sum we want
     * 
     * @return T: The sum of all columns of row r
     */
    T GetRowSum(indextype r);
    
    /**
     * Function to write the matrix content to an ASCII  file
     * 
     *  @param string fname: The name of the file to write
     *  @param char_csep:    The separator character between fields (default: , (comma))
     *  @param withquotes:      bool to indicate if field names in .csv must be written surrounded by quotes.
     * 
     */
    void WriteCsv(std::string fname,char csep=',',bool withquotes=false);
    
    /**
     * Function to write the matrix content to a binary file
     * See format at documentation of JMatrix::WriteBin
     * 
     */
    void WriteBin(std::string fname);
    
    /**
     * Function to get memory in MB used by this symmetric matrix
     * 
     */
    float GetUsedMemoryMB();
    
 private:
     std::vector< std::vector<T> > data;
};

#endif // SYMMETRICMATRIX_H
