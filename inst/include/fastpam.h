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

#ifndef _FASTPAM_H
#define _FASTPAM_H

#include <Rcpp.h>
#include <diftimehelper.h>
#include <symmetricmatrix.h>
#include <threadhelper.h>

// The data types we decide to use for indexes and distances/dissimilarities.
// It is quite unlikely we need to change indextype, unless, to save some memory
// and if you know for sure that there will be no more than 65536 individuals,
// can be defined as unsigned short.
typedef unsigned int indextype;

// The values of these constants are arbitrary, just marks to distinguish them.
// They are meant to mark the different initalization methods
// If you add any other, do at the end and increase the NUM_INIT_METHODS constant.
const unsigned char INIT_METHOD_PREVIOUS=0;
const unsigned char INIT_METHOD_BUILD=1;
const unsigned char INIT_METHOD_LAB=2;
const unsigned char NUM_INIT_METHODS=3;
const std::string init_method_names[NUM_INIT_METHODS]={"PREV","BUILD","LAB"};

// The maximum number of iterations we will allow
const unsigned int  MAX_ITER=1001;

// The maximum number of medoids we allow.
const indextype MAX_MEDOIDS=std::numeric_limits<indextype>::max()-1;

// A convenience constant to indicate a point is not currently assigned to any medoid
const indextype NO_CLUSTER=MAX_MEDOIDS;

// The possible minimum and maximum values the dissimilarity can take. To initialize before loops.
template <typename disttype>
#define MIND std::numeric_limits<disttype>::min()
#define MAXD std::numeric_limits<disttype>::max()
class FastPAM
{
 public:
  FastPAM(SymmetricMatrix<disttype> *Dm,unsigned int num_medoids,unsigned char inimet,int limiter);
  void Init(Rcpp::Nullable<Rcpp::NumericVector> initmedlist,unsigned int nt);
  void Run(unsigned int nt);
  void WriteResultMedoids(std::string medname);
 
  std::vector<indextype> GetMedoids()         { return medoids; };
  std::vector<indextype> GetAssign()          { return nearest; };
  
  std::vector<disttype>  GetTDHistory()       { return TDkeep; };
  std::vector<indextype> GetReassignHistory() { return NpointsChangekeep; };
  
  double GetInTime()  { return(time_in_initialization); };
  double GetOptTime() { return(time_in_optimization); };
  unsigned int GetNumIter() { return(num_iterations_in_opt); };
  
 private:
  SymmetricMatrix<disttype> *D;  // The dissimilarity matrix
  unsigned int nmed;             // The number of medoids we want to find
  indextype num_obs;             // The number of observations; it is equal to the number of rows of D, but just for convenience/clarity.
  unsigned char method;          // The initialization method (see constants to codify methods a few lines up)
  unsigned int maxiter;          // Maximum number of iterations we allow
  bool is_initialized;           // To mark if the chosen initalization algorithm has already be executed.
  
  double time_in_initialization;  // Time in seconds used in the initalization phase (BUILD, ParBUILD or LAB).
  double time_in_optimization;    // Time in seconds used in the optimization phase (FastPAM1 or ParallelFastPAM1).
  unsigned int num_iterations_in_opt;  // Number of itereations used in the optimization phase, never more than maxiter.
  
  // The next fields are filled by initialization (whatever method) and updated by Run
  std::vector<indextype> medoids;     // The current medoids (point index of each one). This is the vector to be returned at the end.
  std::vector<bool>      ismedoid;    // A vector of marks with true for the medoids and false for the others. It is just to accelerate,
                                      // since the medoids vector has already such information.
  std::vector<indextype> nearest;     // The index of the medoid _in_the_array_of_medoids closest to each point
  std::vector<disttype>  dnearest;    // The dissimilarity of every point to its current closest medoid. It plays as a cache
  std::vector<disttype>  dsecond;     // The dissimilarity of every point to its current second closest medoid. It plays as a cache
  
  // These vectors and values are for statistics/information and measures
  disttype               currentTD;          // The value of the optimization function at the current iteration
  std::vector<disttype>  TDkeep;             // Value of TD at each iteration
  indextype              current_npch;       // The value of number of points that have changed cluster at the current iteration
  std::vector<indextype> NpointsChangekeep;  // Number of points that change class at each iteration
 
  // 1) Initialization of the variables; valid for serial and parallel versions
  void InitializeInternals();
  // end 1)
  
  // 2) A function to fill the medoids vector from a previous list; valid for serial and parallel versions
  void InitFromPreviousSet(Rcpp::Nullable<Rcpp::NumericVector> medoidslist);
  // end 2)
   
  // 3) Initalization algorithms
  // 3.1) Brute-force initialization algorithm (BUILD), serial version
  void BUILD();

  // 3.2) Brute-force initialization algorithm, parallel version
  void ParBUILD(unsigned int nt);
       
  // 3.2.1) A structure needed for parallel implementation of BUILD
  struct BUILDThread_IO
  {
      FastPAM *FPp;
      indextype *foundmed;
      disttype *DeltaT;
  };
  // end 3.2.1)
  // 3.2.2) Threads needed for parallel implementation of BUILD
  static void *FindFirstMedoidBUILDThread(void *arg);
  static void *FindSuccessiveMedoidBUILDThread(void *arg);
  // end 3.2.2)
            
  // 3.2) Linear approximative build (LAB), serial version
  void LAB();
  // end 3.2)                  
  // end 3)
  
  // 4) Optimization phase
  // 4.1) Serial version, improved implementation as described in Schubert & Rousseeuw 2021, Algorithm 3
  void RunImprovedFastPAM1();
  // end 4.1)
  
  // 4.2) Parallel version, improved implementation as described in Schubert & Rousseeuw 2021, Algorithm 3
  void RunParallelImprovedFastPAM1(unsigned int nt);
  // 4.2.1) A structure needed for parallel implementation of optimization
  struct FastPAM1Thread_IO
  {
  	FastPAM *FPp;
  	indextype *mst;
  	indextype *xst;
  	indextype *imst;
  	disttype *DeltaTDst;
  	disttype *DeltaTDminusm;         // This will work as an array
  };
  // end 4.2.1)
  // 4.2.2) Thread needed for parallel implementation of optimization
  static void *FastPAM1InternalThread(void *arg);
  // end 4.2.2)
  // end 4.2)
  // end 4)
  
  // 5) Auxiliary functions used inside both versions of optimization
  void FillSecond();
  void SwapRolesAndUpdate(indextype mst,indextype xst,indextype i);
  // end 5)
};


// Prototypes of isolated functions 
std::vector<indextype> randomSample(indextype samplesize, indextype n);
std::vector<indextype> randomSampleExc(indextype samplesize, indextype n,std::vector<bool> &toexclude);
void ApplyPAM(std::string dissim_file,int num_medoids, std::string medoids_file, std::string init_method, std::string PAM_method);
#endif
