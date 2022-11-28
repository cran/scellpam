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

#include <iostream>
#include <sstream>
#include <string>
#include <fullmatrix.h>
#include <sparsematrix.h>
#include <symmetricmatrix.h>
#include <diftimehelper.h>
#include <memhelper.h>
#include <threadhelper.h>
#include <cmath>
#include <unistd.h>

const unsigned char DL1=0x0;   // L1 distance
const unsigned char DL2=0x1;   // L2 distance
const unsigned char DPe=0x2;   // Pearson dissimilarity coefficient

const unsigned char EMPTY=0x0;
const unsigned char IN_FIRST=0x01;
const unsigned char IN_SECOND=0x02;   // This is IN_FIRST << 1
// No need to use this later...
//const unsigned char IN_BOTH=0x03;     // This is IN_FIRST | IN_SECOND

// From now on, and in the .cpp files, counttype is the value type of the input files and disttype the value type of the dissimilarity matrix (our output)
template <typename counttype,typename disttype>
struct args_to_sp_thread
{
    indextype initial_row1;
    indextype final_row1;
    indextype initial_row2;
    indextype final_row2;
    SparseMatrix<counttype> *M;
    SymmetricMatrix<disttype> *D;
    std::vector<counttype> *mu;
    unsigned char dtype;
};

template <typename counttype,typename disttype>
struct args_to_full_thread
{
    unsigned long initial_row1;
    unsigned long final_row1;
    indextype initial_row2;
    indextype final_row2;
    FullMatrix<counttype> *M;
    SymmetricMatrix<disttype> *D;
    std::vector<counttype> *mu;
    unsigned char dtype;
};

// These prototype are declared here so that dissimat.cpp, which calls these functions, can know about them.
template <typename counttype,typename disttype> void CalcAndWriteAuxFull(std::string ifname, std::string ofname, unsigned char dtype, unsigned int nthr,std::string comment="");
template <typename counttype,typename disttype> void CalcAndWriteAuxSparse(std::string ifname, std::string ofname, unsigned char dtype,unsigned int nthr,std::string comment="");

