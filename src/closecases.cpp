/*
 *
 * Copyright (C) 2024 Juan Domingo (Juan.Domingo@uv.es)
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

#include <closecases.h>

extern unsigned char DEB;

// This prototype is declared here so that the call can know it
void CalcAndWriteDissimilarityMatrix(std::string ifname, std::string ofname, std::string distype="L2", std::string restype="float", std::string comment="",int nthreads=0);

void *BasicThreadClose(void *arg)
{
 indextype initial_row = GetField(arg,args_to_close_thread,initial_row);
 indextype final_row = GetField(arg,args_to_close_thread,final_row);
 SymmetricMatrix<float> *M = GetField(arg,args_to_close_thread,M);
 indextype *uq = GetField(arg,args_to_close_thread,uq);
 std::string *method = GetField(arg,args_to_close_thread,method);
 float *value = GetField(arg,args_to_close_thread,value);
 FullMatrix<indextype> *Cl = GetField(arg,args_to_close_thread,Cl);
 
 indextype n=M->GetNRows();
 std::vector<float> thisrow(n);
 std::vector<indextype> idx(n);
 
 for (indextype row=initial_row; row<final_row; row++)
 {
  if ((*method)=="absvalue")
  {
   for (indextype col=0; col<n; col++)
   {
    thisrow[col]=M->Get(row,col);
    // Pearson correlation coefficient is 1-2*pearson_distance. We are interested only in its absolure value.
    thisrow[col]= (thisrow[col]<(*value)) ? 0.0 : fabs(1.0-2*thisrow[col]);
   }
  }
  else
  {
   for (indextype col=0; col<n; col++)
   {
    thisrow[col]=M->Get(row,col);
    // Pearson correlation coefficient is 1-2*pearson_distance. We are interested only in its absolure value.
    thisrow[col]=fabs(1.0-2*thisrow[col]);
   }
  }
  // Here we sort from the highest to the lowest...
  // This is a trick using a vector of indices. iota initializes to 0,1,...,n and the lambda-expression sorts the row and, at the same time, the indices
  iota(idx.begin(), idx.end(), 0);
  std::stable_sort(idx.begin(),idx.end(), 
                   [&thisrow](indextype i1,indextype i2) { return thisrow[i1] > thisrow[i2]; }
                  );
  
  
  // We start storage by 1 since, obviously, the closest to each individual, idx[0], would be itself...
  // Also, indexes are returned as R indexes, so we add 1, and use 0 as mark if the neighbour was discarded becuase its value was below the threshold
  for (indextype e=1; e<(*uq)+1; e++)
   Cl->Set(row,e-1,(thisrow[idx[e]]==0.0) ? 0 : idx[e]+1);
 }
 
 pthread_exit(nullptr);
 // We will never arrive here, but the g++ mingw compiler insists this is a significant warning....
 return nullptr;
}

//' ClosestCases
//'
//' Gets a matrix of size nxp where each file represents an individual (gene, cell,..) and each column its characteristics (counts,...)\cr
//' stored in a file in jmatrix format and returns a matrix of nxp with contains, for each individual the indices of the q individuals\cr
//' closest to it. To do so, the Pearson correlation coefficent is calculated between each row. Then, closeness can be measured directly,\cr
//' as the q values with higher value, but also in a more sophisticated way using the p-values to consider null correlations.
//'
//' @param datafile  A string with the name of the file containing the individuals/characteristics in jmatrix format.
//' @param q         The name of closest related individuals to be returned. Default: 5
//' @param method    The method to be applied to choose the closest values. It must be one of these strings: 'trivial', 'absvalue', 'FDR'.
//' @param dvalue     The value of the absolute value of the Pearson coefficient or of the False Discovery Rate (FDR), depending on the value of the 'method' parameter
//' @param nthreads Number of threads to be used for the parallel calculations with this meaning:\cr
//'                 -1: don't use threads.\cr
//'                  0: let the function choose according to the number of individuals (cells) and to the number of available cores.\cr
//'                  Any possitive number > 1: use that number of threads. You can use even more than cores, but this is discouraged and raises a warning.\cr
//'                 Default: 0.
//' @return          A nxq matrix with the indices (in R notation, starting at 1) of the individuals closest to the individual i at i-th row
//'                  Index will be 0 in some cases if less than q individuals are found to be close according to the 'absvalue' or 'FDR' criteria
//' @examples
//' # To be done
//' @export
// [[Rcpp::export]]
Rcpp::IntegerMatrix ClosestCases(std::string datafile, int q=5, std::string method="trivial", double dvalue=0.0,int nthreads=0)
{
 if (method=="FDR")
 {
  Rcpp::stop("Sorry, FDR method is not yet implemented. We are working on it...\n");
  Rcpp::IntegerMatrix nothing;
  return(nothing);
 }
 if ((method != "trivial") && (method != "absvalue") && (method != "FDR"))
 {
  Rcpp::stop("Parameter method must be one of 'trivial', 'absvalue' or 'FDR'.\n");
  Rcpp::IntegerMatrix nothing;
  return(nothing);
 }

 // Choose a temporaty file. We let this task to R so that the file is erased at the end of the R session
 Rcpp::Environment base("package:base");
 Rcpp::Function GetTmpFile = base["tempfile"];
 Rcpp::CharacterVector rcpp_tdname = GetTmpFile(Rcpp::_["pattern"]="Pdis",Rcpp::_["fileext"]=".bin"); 
 std::string tdname = Rcpp::as<std::string>(rcpp_tdname);

 if (DEB & DEBPP)
  Rcpp::Rcout << "Calculating Pearson dissimilarity matrix, which will be temporary left in " << tdname << "\n";
  
 // We are assuming float precision will be enough for this...
 CalcAndWriteDissimilarityMatrix(datafile,tdname,"Pearson","float","",nthreads);
 
 if (DEB & DEBPP)
  Rcpp::Rcout << "Pearson dissimilarity matrix calculated... ";
 
 // Read the matrix just calculated and stored in file tdname and start playing...
 SymmetricMatrix<float> M(tdname);
 indextype n=M.GetNRows();
 if (DEB & DEBPP)
  Rcpp::Rcout << "and read. It has " << n << " rows.\n";
 
 indextype uq=q;
 FullMatrix<indextype> Cl(n,uq);
 DifftimeHelper Dt;
 
 unsigned int nt=ChooseNumThreads(nthreads);
 if ((n<1000) && (nt!=1))
 {
  nt=1;
  if (DEB & DEBPP)
   Rcpp::Rcout << "We will calculate with a single thread, since you have only " << n << " vectors and the overhead of using threads would be excessive.\n";
 }

 // R passes numeric parameters as double but we are working with floats
 float value=float(dvalue);
 
 // One thread (serial version)
 if (nt==1)
 {
  Dt.StartClock("End of closest neighbour matrix calculation (serial version)."); 
  
  // This will contain each row of the distance matrix, which will have to be ordered
  std::vector<float> thisrow(n);
  std::vector<indextype> idx(n);
 
  for (indextype row=0; row<n; row++)
  {
   if (method=="absvalue")
   {
    for (indextype col=0; col<n; col++)
    {
     thisrow[col]=M.Get(row,col);
     // Pearson correlation coefficient is 1-2*pearson_distance. We are interested only in its absolute value.
     thisrow[col]=(thisrow[col]<value) ? 0.0 : fabs(1.0-2*thisrow[col]);
    }
   }
   else
   {
    for (indextype col=0; col<n; col++)
    {
     thisrow[col]=M.Get(row,col);
     // Pearson correlation coefficient is 1-2*pearson_distance. We are interested only in its absolure value.
     thisrow[col]=fabs(1.0-2*thisrow[col]);
    }
   }
   
   // Here we sort from the highest to the lowest...
   // This is a trick using a vector of indices. iota initializes to 0,1,...,n and the lambda-expression sorts the row and, at the same time, the indices
   iota(idx.begin(), idx.end(), 0);
   std::stable_sort(idx.begin(),idx.end(), 
                    [&thisrow](indextype i1,indextype i2) { return thisrow[i1] > thisrow[i2]; }
                   );
  
   // We start storage by 1 since, obviously, the closest to each individual, idx[0], would be itself...
   // Also, indexes are returned as R indexes, so we add 1, and use 0 as mark if the neighbour was discarded becuase its value was below the threshold
   for (indextype e=1; e<uq+1; e++)
    Cl.Set(row,e-1,(thisrow[idx[e]]==0.0) ? 0 : idx[e]+1);
  }
  Dt.EndClock(DEB & DEBPP);
 }
 else
 {
  Dt.StartClock("End of closest neighbour matrix calculation (parallel version)."); 
  args_to_close_thread *closeargs = new args_to_close_thread [nt];
  
  int rows_per_thread=(n/nt)+1;
  // +1 is because rows are shared "by excess", i.e.: if division is not exact, one row more is assigned
  // to each thread and the last one will have less rows than the others. Just a choice.
  
  for (unsigned int t=0; t<nt; t++)
  {
   closeargs[t].initial_row=t*rows_per_thread;
   closeargs[t].final_row=(t+1)*rows_per_thread;  // The thread is programmed to calculate up to, but not including, the final row
   closeargs[t].M = &M;
   closeargs[t].uq = &uq;
   closeargs[t].method = &method;
   closeargs[t].value = &value;
   closeargs[t].Cl = &Cl;
   if (t==(nt-1))
    closeargs[t].final_row=n;
  }
  
  CreateAndRunThreadsWithDifferentArgs(nt,BasicThreadClose,(void *)closeargs,sizeof(args_to_close_thread));
  
  delete[] closeargs;
  Dt.EndClock(DEB & DEBPP);
 }
 
 Rcpp::IntegerMatrix ret(n,uq);
 for (indextype row=0; row<n; row++)
  for (indextype col=0; col<uq; col++)
   ret(row,col)=Cl.Get(row,col);
   
 return(ret);
}

