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

#include <Rcpp.h>
#include <thread>
#include <diftimehelper.h>
#include <fullmatrix.h>
#include <sparsematrix.h>
#include <symmetricmatrix.h>
#include <threadhelper.h>
#include <memhelper.h>
#include <debugpar.h>

// The values of the silhouette will be stored as double. Since its number is always linear
// with the number of points, this should not increase memoary usage too much and it is simpler.
typedef double siltype;

// A vector of structures like the following one will be created, one structure per point
// Its objective is to be able to return to R a matrix compatible with package cluster which allow representation in standard form.
typedef struct
    {
     indextype pnum;                    // The original point number (to keep it after sorting)
     indextype ownclus;			// The cluster the point belongs to (C++-numbering, 0..nclus-1)
     indextype neiclus;                 // The cluster which is at minimal average distance (except the own cluster), i.e.: the closest neighbour (C++-numbering)
     siltype  silvalue;                 // The silhouette value
    } silinfo;
    
template <typename disttype>
struct SilhoutteThread_IO
    {
     indextype num_obs;
     indextype nmed;
     const std::vector<indextype> *nearest;
     std::vector<siltype> *current_sil;
     std::vector<unsigned long> *hist;
     std::vector<silinfo> *silres;
     SymmetricMatrix<disttype> *D;
    };

