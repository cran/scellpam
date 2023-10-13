/*
 *
 * Copyright (C) 2023 Juan Domingo (Juan.Domingo@uv.es)
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

#include <dissimmat.h>

extern unsigned char DEB;

// This function will fill part of the distance matrix D, concretely, lines between initial_row and (but not including) final_row
template <typename counttype,typename disttype>
void FillMetricMatrixFromSparse(indextype initial_row,indextype final_row,SparseMatrix<counttype> *M,SymmetricMatrix<disttype> *D,bool L1dist)
{
 disttype d,dif;
 indextype ncols=M->GetNCols();
 indextype nrows=D->GetNRows();
 
 // This should not be done from inside a thread. But, if we have failed anyway....
 if ( (initial_row >= nrows) || (final_row > nrows) )
 {
     std::ostringstream errst;
     errst << "Error in FillMetricMatrixFromSparse: either start of area at " << initial_row << " or end of area at " << final_row << " or both are outside matrix limits.\n";
     Rcpp::stop(errst.str());
     return;
 }
 
 counttype *va = new counttype [ncols];
 counttype *vb = new counttype [ncols];
 unsigned char *mark = new unsigned char [ncols];
 unsigned char *rowmark = new unsigned char [ncols];
 
 for (indextype rowA=initial_row; rowA<final_row; rowA++)
 {
  // The values of the current row (let's call it rowA) and the places they are, are stored in va and rowmark respectively
  memset((void *)va,0x0,ncols*sizeof(counttype));
  memset((void *)rowmark,EMPTY,ncols);
  
  M->GetSparseRow(rowA,rowmark,IN_FIRST,va);
  
  // The next loop calculates the distance between the current row (rowA) and all others with numbers below its own number, let's call rowB to each.
  // (remember that to fill distance matrix we only need to fill the lower-diagonal part...)
  
  for (unsigned rowB=0; rowB<rowA; rowB++)
  {
   // The mark array is initialized to the marks of rowA
   memcpy(mark,rowmark,ncols);
   
   // and vb will contains the values of rowB. First, we clean it...
   memset((void *)vb,0x0,ncols*sizeof(counttype));
   
   // and now, we fill its values and update mark adding 2 to the rigth places. In that way, mark will return with 0, 1, 2 or 3
   // for the cases of "No value at that place in any row", "Value only in first row", "Value only in second row" or "Value in both", respectively.
   M->GetSparseRow(rowB,mark,IN_SECOND,vb);
   
   d=0.0;
   for (unsigned col=0; col<ncols; col++)
   {
    // That component of the vector is updated only if there is something to consider...
    if (mark[col]!=EMPTY)
    {
        dif = ((mark[col]==IN_FIRST) ? disttype(va[col]) : ((mark[col]==IN_SECOND) ? -disttype(vb[col]) : disttype(va[col])-disttype(vb[col]) ));
        d += (L1dist ? fabs(dif) : dif*dif);
    }
   } 
   
   D->Set(rowA,rowB,(L1dist ? d : sqrt(d)));
  }
  
  // This is just to set the main diagonal.
  D->Set(rowA,rowA,disttype(0));
 }
 
 delete[] va;
 delete[] vb;
 delete[] mark;
 delete[] rowmark;
}

template void FillMetricMatrixFromSparse(indextype initial_row,indextype final_row,SparseMatrix<float> *M,SymmetricMatrix<float> *D,bool L1dist);
template void FillMetricMatrixFromSparse(indextype initial_row,indextype final_row,SparseMatrix<double> *M,SymmetricMatrix<float> *D,bool L1dist);
template void FillMetricMatrixFromSparse(indextype initial_row,indextype final_row,SparseMatrix<float> *M,SymmetricMatrix<double> *D,bool L1dist);
template void FillMetricMatrixFromSparse(indextype initial_row,indextype final_row,SparseMatrix<double> *M,SymmetricMatrix<double> *D,bool L1dist);

// This function will fill part of the distance matrix D, concretely, lines between initial_row and (but not including) final_row
template <typename counttype,typename disttype>
void FillCosMatrixFromSparse(indextype initial_row,indextype final_row,SparseMatrix<counttype> *M,SymmetricMatrix<disttype> *D)
{
 disttype sc1,sc2,prod;
 indextype ncols=M->GetNCols();
 indextype nrows=D->GetNRows();
 
 // This should not be done from inside a thread. But, if we have failed anyway....
 if ( (initial_row >= nrows) || (final_row > nrows) )
 {
     std::ostringstream errst;
     errst << "Error in FillCosMatrixFromFull: either start of area at " << initial_row << " or end of area at " << final_row << " or both are outside matrix limits.\n";
     Rcpp::stop(errst.str());
     return;
 }
 
 counttype *va = new counttype [ncols];
 counttype *vb = new counttype [ncols];
 unsigned char *mark = new unsigned char [ncols];
 unsigned char *rowmark = new unsigned char [ncols];
 
 for (indextype rowA=initial_row; rowA<final_row; rowA++)
 {
  // The values of the current row (let's call it rowA) and the places they are, are stored in va and rowmark respectively
  memset((void *)va,0x0,ncols*sizeof(counttype));
  memset((void *)rowmark,EMPTY,ncols);
  
  M->GetSparseRow(rowA,rowmark,IN_FIRST,va);
  
  // The next loop calculates the distance between the current row (rowA) and all others with numbers below its own number, let's call rowB to each.
  // (remember that to fill distance matrix we only need to fill the lower-diagonal part...)
  
  for (unsigned rowB=0; rowB<rowA; rowB++)
  {
   // The mark array is initialized to the marks of rowA
   memcpy(mark,rowmark,ncols);
   
   // and vb will contain the values of rowB. First, we clean it...
   memset((void *)vb,0x0,ncols*sizeof(counttype));
   
   // and now, we fill its values and update mark adding 2 to the right places. In that way, mark will return with 0, 1, 2 or 3
   // for the cases of "No value at that place in any row", "Value only in first row", "Value only in second row" or "Value in both", respectively.
   M->GetSparseRow(rowB,mark,IN_SECOND,vb);
   
   sc1=sc2=prod=disttype(0);
   for (unsigned col=0; col<ncols; col++)
   {
    if (mark[col])
    {
     switch (mark[col])
     {
      case IN_FIRST: sc1 += disttype(va[col]*va[col]);
                     break;
      case IN_SECOND: sc2 += disttype(vb[col]*vb[col]);
                      break;
      case IN_BOTH: sc1 += disttype(va[col]*va[col]);
                    sc2 += disttype(vb[col]*vb[col]); 
                    prod += disttype(va[col]*vb[col]);
                    break;
      default: break;
     }
    }              
   } 
   // Precision problems in the sqrt calculation might yield an extremely small negative result. we clip it to zero
   disttype val=1.0-(prod/(sqrt(sc1)*sqrt(sc2)));
   D->Set(rowA,rowB,val<0.0 ? 0.0 : val);
  }
  
  // This is just to set the main diagonal.
  D->Set(rowA,rowA,disttype(0));
 }
 
 delete[] va;
 delete[] vb;
 delete[] mark;
 delete[] rowmark;
}

template void FillCosMatrixFromSparse(indextype initial_row,indextype final_row,SparseMatrix<float> *M,SymmetricMatrix<float> *D);
template void FillCosMatrixFromSparse(indextype initial_row,indextype final_row,SparseMatrix<double> *M,SymmetricMatrix<float> *D);
template void FillCosMatrixFromSparse(indextype initial_row,indextype final_row,SparseMatrix<float> *M,SymmetricMatrix<double> *D);
template void FillCosMatrixFromSparse(indextype initial_row,indextype final_row,SparseMatrix<double> *M,SymmetricMatrix<double> *D);
       
template <typename counttype,typename disttype>                    
void CalculateMeansFromSparse(SparseMatrix<counttype> &M,std::vector<disttype> &mu)
{
 indextype ncells=M.GetNRows();
 indextype ngenes=M.GetNCols();
 
 disttype s;   
 for (indextype gene=0; gene<ngenes; gene++)
 {
  s=disttype(0); 
  for (indextype cell=0; cell<ncells; cell++)
   s += disttype(M.Get(cell,gene));
  mu.push_back(s/disttype(ncells)); 
 }
}

template void CalculateMeansFromSparse(SparseMatrix<float> &M,std::vector<float> &mu);
template void CalculateMeansFromSparse(SparseMatrix<float> &M,std::vector<double> &mu);
template void CalculateMeansFromSparse(SparseMatrix<double> &M,std::vector<float> &mu);
template void CalculateMeansFromSparse(SparseMatrix<double> &M,std::vector<double> &mu);

template <typename counttype,typename disttype>
void CalculateVariancesFromSparse(SparseMatrix<counttype> &M,std::vector<disttype> &mu,std::vector<disttype> &cvar)
{
 indextype ncells=M.GetNRows();
 indextype ngenes=M.GetNCols();
 
 disttype st,dif;
 for (indextype gene=0; gene<ngenes; gene++)
 {
  st=disttype(0); 
  for (indextype cell=0; cell<ncells; cell++)
  {
   dif = disttype(M.Get(cell,gene))-disttype(mu[gene]);
   st += dif*dif;
  }
  cvar.push_back(st/disttype(ncells-1)); 
 }
}

template void CalculateVariancesFromSparse(SparseMatrix<float> &M,std::vector<float> &mu,std::vector<float> &cvar);
template void CalculateVariancesFromSparse(SparseMatrix<double> &M,std::vector<float> &mu,std::vector<float> &cvar);
template void CalculateVariancesFromSparse(SparseMatrix<float> &M,std::vector<double> &mu,std::vector<double> &cvar);
template void CalculateVariancesFromSparse(SparseMatrix<double> &M,std::vector<double> &mu,std::vector<double> &cvar);

// This function will fill part of the distance matrix D, concretely, lines between initial_row and (but not including) final_row
template <typename counttype,typename disttype>
void FillPearsonMatrixFromSparse(indextype initial_row,indextype final_row,SparseMatrix<counttype> *M,std::vector<disttype> *mu,SymmetricMatrix<disttype> *D)
{
 disttype da,db,sxx,syy,sxy,den,pearson;
 disttype dtol=std::numeric_limits<disttype>::epsilon();
 indextype ncols=M->GetNCols();
 indextype nrows=D->GetNRows();
 
 // This should not be done from inside a thread. But, if we have failed anyway....
 if ( (initial_row >= nrows) || (final_row > nrows) )
 {
     std::ostringstream errst;
     errst << "Error in FillPearsonMatrixFromSparse: either start of area at " << initial_row << " or end of area at " << final_row << " or both are outside matrix limits.\n";
     Rcpp::stop(errst.str());
     return;
 }
 
 // Memory is booked for all the row. Only part of it (from initial_col to initial_col+extension) will be
 // used, but we have memory enough and it is easier this way. The same applies to the arrays of marks.
 counttype *va = new counttype [ncols];
 counttype *vb = new counttype [ncols];
 
 for (unsigned rowA=initial_row; rowA<final_row; rowA++)
 {
  // The values of the current row (let's call it rowA) are stored in va 
  memset((void *)va,0x0,ncols*sizeof(counttype));
  
  M->GetRow(rowA,va);
  
  // The next loop calculates the distance between the current row (rowA) and all others with numbers below its own number, let's call rowB to each.
  // (remember that to fill distance matrix we only need to fill the lower-diagonal part...)
  
  for (unsigned rowB=0; rowB<rowA; rowB++)
  {
   // bb will contains the values of rowB. First, we clean it...
   memset((void *)vb,0x0,ncols*sizeof(counttype));
   
   // and now, we fill its values
   M->GetRow(rowB,vb);
    
   sxx=sxy=syy=0.0;
   for (unsigned col=0; col<ncols; col++)
   {
    da=disttype(va[col])-((*mu)[col]);
    db=disttype(vb[col])-((*mu)[col]);
    sxx += (da*da);
    syy += (db*db);
    sxy += (da*db);
   }
   
   // This is to prevent from numerical problems
   den=disttype(sqrt(double(sxx)))*disttype(sqrt(double(syy)));
   if (den==0.0)   // This is the pathological case in which both vectors are the null vector. Then, they are of course completely "similar" (indeed, identical...)
    D->Set(rowA,rowB,disttype(0.0));
   else
   {
    pearson=0.5-(sxy/den/2.0);
    D->Set(rowA,rowB,(fabs(pearson)<dtol) ? disttype(0.0) : pearson);
   }
  }
  // This is just to set the main diagonal.
  D->Set(rowA,rowA,disttype(0));
 }
 delete[] va;
 delete[] vb;
}

template void FillPearsonMatrixFromSparse(indextype initial_row,indextype final_row,SparseMatrix<float> *M,std::vector<float> *mu,SymmetricMatrix<float> *D);
template void FillPearsonMatrixFromSparse(indextype initial_row,indextype final_row,SparseMatrix<double> *M,std::vector<float> *mu,SymmetricMatrix<float> *D);
template void FillPearsonMatrixFromSparse(indextype initial_row,indextype final_row,SparseMatrix<float> *M,std::vector<double> *mu,SymmetricMatrix<double> *D);
template void FillPearsonMatrixFromSparse(indextype initial_row,indextype final_row,SparseMatrix<double> *M,std::vector<double> *mu,SymmetricMatrix<double> *D);

// This function will fill part of the distance matrix D, concretely, lines between initial_row and (but not including) final_row
template <typename counttype,typename disttype>
void FillWEucMatrixFromSparse(indextype initial_row,indextype final_row,SparseMatrix<counttype> *M,std::vector<disttype> *cvar,SymmetricMatrix<disttype> *D)
{
 disttype d,dif;
 indextype ncols=M->GetNCols();
 indextype nrows=D->GetNRows();
 
 // This should not be done from inside a thread. But, if we have failed anyway....
 if ( (initial_row >= nrows) || (final_row > nrows) )
 {
     std::ostringstream errst;
     errst << "Error in FillWEucMatrixFromSparse: either start of area at " << initial_row << " or end of area at " << final_row << " or both are outside matrix limits.\n";
     Rcpp::stop(errst.str());
     return;
 }
 
 counttype *va = new counttype [ncols];
 counttype *vb = new counttype [ncols];
 unsigned char *mark = new unsigned char [ncols];
 unsigned char *rowmark = new unsigned char [ncols];

 for (indextype rowA=initial_row; rowA<final_row; rowA++)
 {
  // The values of the current row (let's call it rowA) and the places they are, are stored in va and rowmark respectively
  memset((void *)va,0x0,ncols*sizeof(counttype));
  memset((void *)rowmark,EMPTY,ncols);
  
  M->GetSparseRow(rowA,rowmark,IN_FIRST,va);
  
  // The next loop calculates the distance between the current row (rowA) and all others with numbers below its own number, let's call rowB to each.
  // (remember that to fill distance matrix we only need to fill the lower-diagonal part...)
  
  for (indextype rowB=0; rowB<rowA; rowB++)
  {
   // The mark array is initialized to the marks of rowA
   memcpy(mark,rowmark,ncols);
   
   // and vb will contain the values of rowB. First, we clean it...
   memset((void *)vb,0x0,ncols*sizeof(counttype));
   
   // and now, we fill its values and update mark adding 2 to the rigth places. In that way, mark will return with 0, 1, 2 or 3
   // for the cases of "No value at that place in any row", "Value only in first row", "Value only in second row" or "Value in both", respectively.
   M->GetSparseRow(rowB,mark,IN_SECOND,vb);
   
   d=0.0;
   for (unsigned col=0; col<ncols; col++)
   {
    // That component of the vector is updated only if there is something to consider...
    if (mark[col]!=EMPTY)
    {
        dif = ((mark[col]==IN_FIRST) ? disttype(va[col]) : ((mark[col]==IN_SECOND) ? -disttype(vb[col]) : disttype(va[col])-disttype(vb[col]) ));
        // The case of null variance does not happen: it would mean all the column is null, but in that case mark[col] would be EMPTY and we wouldn't be inside this if
        d += (dif*dif/(*cvar)[col]);
    }
   } 
   
   D->Set(rowA,rowB,sqrt(d));
  }
  
  // This is just to set the main diagonal.
  D->Set(rowA,rowA,disttype(0));
 }
 delete[] va;
 delete[] vb;
 delete[] mark;
 delete[] rowmark;
}

template void FillWEucMatrixFromSparse(indextype initial_row,indextype final_row,SparseMatrix<float> *M,std::vector<float> *cvar,SymmetricMatrix<float> *D);
template void FillWEucMatrixFromSparse(indextype initial_row,indextype final_row,SparseMatrix<double> *M,std::vector<float> *cvar,SymmetricMatrix<float> *D);
template void FillWEucMatrixFromSparse(indextype initial_row,indextype final_row,SparseMatrix<float> *M,std::vector<double> *cvar,SymmetricMatrix<double> *D);
template void FillWEucMatrixFromSparse(indextype initial_row,indextype final_row,SparseMatrix<double> *M,std::vector<double> *cvar,SymmetricMatrix<double> *D);

template <typename counttype,typename disttype>
void *BasicThreadSparse(void *arg)
{
 indextype initial_row1 = GetFieldDT(arg,args_to_sp_thread,counttype,disttype,initial_row1);
 indextype final_row1 = GetFieldDT(arg,args_to_sp_thread,counttype,disttype,final_row1);
 indextype initial_row2 = GetFieldDT(arg,args_to_sp_thread,counttype,disttype,initial_row2);
 indextype final_row2 = GetFieldDT(arg,args_to_sp_thread,counttype,disttype,final_row2);
 SparseMatrix<counttype> *M = GetFieldDT(arg,args_to_sp_thread,counttype,disttype,M);
 SymmetricMatrix<disttype> *D = GetFieldDT(arg,args_to_sp_thread,counttype,disttype,D);
 std::vector<disttype> *mu_or_var = GetFieldDT(arg,args_to_sp_thread,counttype,disttype,muvar);
 unsigned char dtype = GetFieldDT(arg,args_to_sp_thread,counttype,disttype,dtype);

 switch (dtype)
 {
  case DL1: {
  		FillMetricMatrixFromSparse(initial_row1,final_row1,M,D,true);
            	FillMetricMatrixFromSparse(initial_row2,final_row2,M,D,true);
            }
            break;
  case DL2: {
  		FillMetricMatrixFromSparse(initial_row1,final_row1,M,D,false);
            	FillMetricMatrixFromSparse(initial_row2,final_row2,M,D,false);
            }
            break;
  case DPe: {
  		FillPearsonMatrixFromSparse(initial_row1,final_row1,M,mu_or_var,D);
            	FillPearsonMatrixFromSparse(initial_row2,final_row2,M,mu_or_var,D);
            }
            break;
  case DCo: {
  		FillCosMatrixFromSparse(initial_row1,final_row1,M,D);
            	FillCosMatrixFromSparse(initial_row2,final_row2,M,D);
            }
            break;
  case DWe: {
  		FillWEucMatrixFromSparse(initial_row1,final_row1,M,mu_or_var,D);
            	FillWEucMatrixFromSparse(initial_row2,final_row2,M,mu_or_var,D);
            }
            break;
  default: break;
 }
 pthread_exit(nullptr);
 // We will never arrive here, but the g++ mingw compiler insists this is a significant warning....
 return nullptr;   
}

template void *BasicThreadSparse<float,float>(void *arg);
template void *BasicThreadSparse<float,double>(void *arg);
template void *BasicThreadSparse<double,float>(void *arg);
template void *BasicThreadSparse<double,double>(void *arg);

template <typename counttype,typename disttype>
void CalcAndWriteAuxSparse(std::string ifname, std::string ofname, unsigned char dtype,unsigned int nthr,std::string comment)
{
 SparseMatrix<counttype> M(ifname);
 indextype nrows=M.GetNRows();
 if (DEB & DEBPP)
 {
  Rcpp::Rcout << "Read sparse matrix from file " << ifname << ".\n   Its size is [" << M.GetNRows() << " x " << M.GetNCols() << "] and it uses " << M.GetUsedMemoryMB() << " MBytes.\n";
  Rcpp::Rcout << "Creating dissimilarity matrix of size (" << nrows << "x" << nrows << ")\n";
 }
 
 MemoryWarnings(nrows,sizeof(disttype));
  
 SymmetricMatrix<disttype> D(nrows); 
 
 if (DEB & DEBPP)
 {
  std::string dummy;
  switch (dtype)
  {
   case DL1: dummy="L1"; break;
   case DL2: dummy="L2"; break;
   case DPe: dummy="Pearson"; break;
   case DCo: dummy="Cosine"; break;
   case DWe: dummy="WeightedEuclidean"; break;
  }
 }
 
 std::vector<disttype> mu,cvar;
 if ((dtype==DPe) || (dtype==DWe))
 {
  if (DEB & DEBPP)
   Rcpp::Rcout << "Calculating vector of means used by Pearson dissimilarity and weighted Euclidean distance...\n";
  CalculateMeansFromSparse(M,mu);
  if (mu.size()!=M.GetNCols())
   Rcpp::stop("Error from CalcAndWriteAuxSparse: length of vector of means is not the number of columns of the data matrix.\n");
 }
 if (dtype==DWe)
 {
  if (DEB & DEBPP)
   Rcpp::Rcout << "Calculating vector of variances used by weigthed Euclidean distance...\n";
  CalculateVariancesFromSparse(M,mu,cvar);
  if (cvar.size()!=M.GetNCols())
   Rcpp::stop("Error from CalcAndWriteAuxFull: length of vector of variances is not the number of columns of the data matrix.\n");
 }
 
 if ((nrows<1000) && (nthr!=1))
 {
  nthr=1;
  if (DEB & DEBPP)
   Rcpp::Rcout << "Calculating with a single thread, since you have only " << nrows << " vectors and the overhead of using threads would be excessive.\n";
 }
 
 DifftimeHelper Dt;
 
 if (nthr==1)
 {
  Dt.StartClock("End of dissimilarity matrix calculation (serial version).");
  switch (dtype)
  {
   case DL1: FillMetricMatrixFromSparse(0,D.GetNRows(),&M,&D,true); break;
   case DL2: FillMetricMatrixFromSparse(0,D.GetNRows(),&M,&D,false); break;
   case DPe: FillPearsonMatrixFromSparse(0,D.GetNRows(),&M,&mu,&D); break;
   case DCo: FillCosMatrixFromSparse(0,D.GetNRows(),&M,&D); break;
   case DWe: FillWEucMatrixFromSparse(0,D.GetNRows(),&M,&cvar,&D); break;
   default: break;
  }
  Dt.EndClock(DEB & DEBPP);
 }
 else
 {
  Dt.StartClock("End of dissimilarity matrix calculation (parallel version)."); 

  indextype nrows=D.GetNRows();
  
  args_to_sp_thread<counttype,disttype> *spargs = new args_to_sp_thread<counttype,disttype> [nthr];
  
  // This strange distribution of rows among threads (each thread has two intervals which are symmetric with respect to the middle of the matrix rows)
  // is chosen to balance as much as possible the number of distances calculated by each thread. For each 'short' row (those before the middle) the same
  // thread does also a 'long' row (those after the middle)
  for (unsigned int t=0; t<nthr; t++)
  {
   spargs[t].initial_row1 = indextype(float(t)*float(nrows)/(2*float(nthr)));
   if (t==nthr-1)
    spargs[t].final_row1 = nrows/2;
   else
    spargs[t].final_row1 = indextype(float(t+1)*float(nrows)/(2*float(nthr)));
      
   spargs[t].initial_row2 = nrows-spargs[t].final_row1 - 1;
   if ((t==nthr-1) && (!(nrows%2)))
    spargs[t].initial_row2++;
    
   spargs[t].final_row2 = nrows-spargs[t].initial_row1 - 1;
   if (t==0)
    spargs[t].final_row2++;
    
   spargs[t].M = &M;
   spargs[t].D = &D;
   switch (dtype)
   {
    case DPe: spargs[t].muvar = &mu; break;
    case DWe: spargs[t].muvar = &cvar;  break;
    case DL1: spargs[t].muvar = nullptr; break;
    case DL2: spargs[t].muvar = nullptr; break;
    case DCo: spargs[t].muvar = nullptr; break;
    default: break;
   }
   spargs[t].dtype = dtype;
  }
  if (DEB & DEBPP)
  {
     Rcpp::Rcout << "Launching up to " << nthr << " simultaneous threads.\n";
     /*
     Rcpp::Rcout << "Starting/end rows for each thread:\n";
     for (unsigned int t=0; t<nthr; t++)
      Rcpp::Rcout << t << ": [" << spargs[t].initial_row1 << "," << spargs[t].final_row1 << "]+[" << spargs[t].initial_row2 << "," << spargs[t].final_row2 << "]; ";
     Rcpp::Rcout << "\n";
     */
  }
  
  CreateAndRunThreadsWithDifferentArgs(nthr,BasicThreadSparse<counttype,disttype>,(void *)spargs,sizeof(args_to_sp_thread<counttype,disttype>));
  
  delete[] spargs;
 
  Dt.EndClock(DEB & DEBPP);
 }
 
 std::vector<std::string> rnames=M.GetRowNames();
 if (!rnames.empty())     
  D.SetRowNames(rnames);
  
 if (comment!="")
  D.SetComment(comment);
     

 D.WriteBin(ofname);
 if (DEB & DEBPP)
     Rcpp::Rcout << "Output binary file " << ofname << " written.\n";
}

template void CalcAndWriteAuxSparse<float,float>(std::string ifname, std::string ofname, unsigned char dtype,unsigned int nthr,std::string comment="");
template void CalcAndWriteAuxSparse<float,double>(std::string ifname, std::string ofname, unsigned char dtype,unsigned int nthr,std::string comment="");
template void CalcAndWriteAuxSparse<double,float>(std::string ifname, std::string ofname, unsigned char dtype,unsigned int nthr,std::string comment="");
template void CalcAndWriteAuxSparse<double,double>(std::string ifname, std::string ofname, unsigned char dtype,unsigned int nthr,std::string comment="");


