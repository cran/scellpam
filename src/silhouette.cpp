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

#include<silhouette.h>
#include<pthread.h>

extern unsigned char DEB;

// Auxiliary function to sort the silhouette as needed to make an object suitable for the cluster package
void SortByClusterAndSilhouette(std::vector<silinfo> &silres)
{
 sort(silres.begin(),silres.end(),[](silinfo a,silinfo b)
 {
  if (a.ownclus < b.ownclus)
   return(true);
  else
   if (a.ownclus==b.ownclus)
    return(a.silvalue > b.silvalue);
   else
    return(false);
 });
}    

// Thread used to calculate the silhouette in its parallel implementation
template <typename disttype>
void *SilhoutteThread(void *arg)
{
  unsigned int numthreads = GetNumThreads(arg);
  unsigned int current_thread_num = GetThisThreadNumber(arg);
  
  indextype num_obs = GetField(arg,SilhoutteThread_IO<disttype>,num_obs);
  indextype nmed = GetField(arg,SilhoutteThread_IO<disttype>,nmed);
  const std::vector<indextype> *nearest = GetField(arg,SilhoutteThread_IO<disttype>,nearest);
  std::vector<siltype> *current_sil = GetField(arg,SilhoutteThread_IO<disttype>,current_sil);
  std::vector<unsigned long> *hist = GetField(arg,SilhoutteThread_IO<disttype>,hist);
  std::vector<silinfo> *silres = GetField(arg,SilhoutteThread_IO<disttype>,silres);
  SymmetricMatrix<disttype> *D = GetField(arg,SilhoutteThread_IO<disttype>,D);
  
  indextype start = current_thread_num*(num_obs/numthreads);
  indextype end   = (current_thread_num == numthreads-1) ? num_obs : (current_thread_num+1)*(num_obs/numthreads);
  
  // To interpret all variables and operation here, please look at the comments in the serial version below
  siltype *a = new siltype [end-start];
  siltype *b = new siltype [end-start];
  siltype *bav = new siltype [nmed];
  siltype dmin;
  indextype which_neimin=nmed+1;
  for (indextype q=start; q<end; q++)
  {
     if ((*hist)[(*nearest)[q]]==1)
      (*current_sil)[q]=0.0;
     else
     {
      for (indextype m=0; m<nmed; m++)
       bav[m]=0.0;
      
      for (indextype q1=0; q1<num_obs; q1++)
       bav[(*nearest)[q1]] += D->Get(q,q1);
          
      for (indextype m=0; m<nmed; m++)
       if (m==(*nearest)[q])
        bav[m] /= double((*hist)[m]-1);
       else
        bav[m] /= double((*hist)[m]);
      
      a[q-start] = bav[(*nearest)[q]];
        
      dmin=std::numeric_limits<siltype>::max();
      for (indextype m=0; m<nmed; m++)
       if ( (m!=(*nearest)[q]) && (bav[m]<dmin) )
       {
           which_neimin=m;
           dmin = bav[m];
       } 
       
      b[q-start]=dmin;   
      
      (*current_sil)[q]=(b[q-start]-a[q-start])/std::max(a[q-start],b[q-start]);
     }
     
     (*silres)[q].neiclus=which_neimin;
     (*silres)[q].silvalue=(*current_sil)[q];
  }
 
  delete[] a;
  delete[] b;
  delete[] bav;
  
  pthread_exit(nullptr);
  // We will never arrive here, but the g++ mingw compiler insists this is a significant warning....
  return nullptr;  
}

template void *SilhoutteThread<float>(void *arg);
template void *SilhoutteThread<double>(void *arg);

// Auxiliary function with the real implementation of silhouette (serial version)
template <typename disttype>
void SilhouetteSerial(indextype num_obs,indextype nmed,
                      std::vector<indextype> &nearest,std::vector<siltype> &current_sil,std::vector<unsigned long> &hist,
                      std::vector<silinfo> &silres,SymmetricMatrix<disttype> &D)
{
 siltype *a = new siltype [num_obs];
 siltype *b = new siltype [num_obs];
 siltype *bav = new siltype [nmed];
 siltype dmin;
    
 indextype which_neimin=nmed+1;
 for (indextype q=0; q<num_obs; q++)
 {
   // Special case: cluster with one isolated point. Silhouette in this case is defined as 0
   if (hist[nearest[q]]==1)
     current_sil[q]=0.0;
   else
   {  
    // bav will contain the average distance between this points and the others in each cluster, including its own cluster
    for (indextype m=0; m<nmed; m++)
     bav[m]=0.0;
         
    for (indextype q1=0; q1<num_obs; q1++)
     bav[nearest[q1]] += D.Get(q,q1);
          
    // bav contains now the sum of distances to point q, by cluster. Let's divide to calculate the average.
    // The 'minus 1' is because the point itself is not counted. This is the definition of silhouette.
    // In this case, the denominator cannot be 0, since the special case of hist[nearest[q]]==1 was managed before
    // Also, it is not possible to have hist==0 for other clusters, since every cluster is requested to have at least one point
    for (indextype m=0; m<nmed; m++)
     if (m==nearest[q])
      bav[m] /= double(hist[m]-1);
     else
      bav[m] /= double(hist[m]);
           
    // Now, for each point q, a will contain the average distance to the other points in its own cluster:
    a[q] = bav[nearest[q]];
         
    // and b will contain the minimal average distance to points in _other_ clusters.
    // This is why we leave out the own cluster of the point
    // We keep also the number of the cluster at minimal distance, in which_neimin
    dmin=std::numeric_limits<siltype>::max();
    for (indextype m=0; m<nmed; m++)
     if ( (m!=nearest[q]) && (bav[m]<dmin) )
     {
      which_neimin=m;
      dmin = bav[m];
     }  
     
    b[q]=dmin;   
          
    current_sil[q]=(b[q]-a[q])/std::max(a[q],b[q]);
   }
        
   silres[q].neiclus=which_neimin;
   silres[q].silvalue=current_sil[q];       
 }
 
 delete[] a;
 delete[] b;
 delete[] bav;
  
}

// Auxiliary function to calculate silhouette. It implements both versions inside: serial an parallel
template <typename disttype>
Rcpp::NumericVector CalculateSilhouetteAux(Rcpp::NumericVector cl,std::string fdist,unsigned int nt)
{
 DifftimeHelper Dt;
 if (nt==1)
 {
   if (DEB & DEBPP)
   {
    Rcpp::Rcout << "   Calculating silhouette (serial implementation)...\n";
    Rcpp::Rcout.flush();
   }
   Dt.StartClock("Finished serial implementation of silhouette (including dissimilarity matrix load).");
 }
 else
 {
   if (DEB & DEBPP)
   {
    Rcpp::Rcout << "   Calculating silhouette (parallel version) with " << nt << " threads.\n";
    Rcpp::Rcout.flush();
   }
   Dt.StartClock("Finished parallel implementation of silhouette (including dissimilarity matrix load).");
 }
    
 SymmetricMatrix<disttype> D(fdist);

 indextype num_obs=D.GetNRows();
 
 // Nearest is the vector of nearest cluster number (in C++ numbering, from 0)
 std::vector<indextype> nearest;

 indextype mincl=num_obs;
 indextype maxcl=0;
 
 // This is just to check cl vector, which is in R numbering (from 1) and set nearest from it
 for (size_t t=0; t<size_t(cl.length()); t++)
 {
  if (cl[t]<1 || cl[t]>num_obs)
   Rcpp::stop("The clasification array contains at least one invalid value (outside the range 1.. number_of_points).\n");
  if (cl[t]-1 < mincl)
   mincl=cl[t]-1;
  if (cl[t]-1 > maxcl)
   maxcl=cl[t]-1;
  nearest.push_back(cl[t]-1);
 }
 if (mincl!=0)
  Rcpp::stop("The classification array has not 1 as minimum value.\n");
  
 // The number of clusters (number of medoids)
 indextype nmed=maxcl+1;
 
 if (num_obs!=nearest.size())
  Rcpp::stop("Different number of points in the array of classes and in the dissimilarity matrix.\n");
 
 if (DEB & DEBPP)
  Rcpp::Rcout << num_obs << " points classified in " << nmed << " classes.\n";
 
 silinfo sdummy;
 std::vector<silinfo> silres;
 
 // hist is the histogram to know how many individuals are assigned to each cluster
 std::vector<unsigned long> hist;
 hist.resize(nmed,0);
 for (indextype q=0; q<num_obs; q++)
 {
  hist[nearest[q]]++;
  // The structure with the silhouette info for each point is initialized: original point number and its cluster number, which will not be changed.
  sdummy.pnum=q;
  sdummy.ownclus=nearest[q];
  // Also, the nearest (other than its own) cluster and silhouette value are initalized to control values, to be filled.
  sdummy.neiclus=nmed;            // Absurd values to serve as control..
  sdummy.silvalue=std::numeric_limits<siltype>::max();
  silres.push_back(sdummy);
 }
  
 std::vector<siltype> current_sil;
 current_sil.resize(num_obs,siltype(0));
 
 if (nt==1)
    SilhouetteSerial(num_obs,nmed,nearest,current_sil,hist,silres,D);   
 else
 {
    SilhoutteThread_IO<disttype> *silargs = new SilhoutteThread_IO<disttype> [nt];
    for (unsigned int t=0; t<nt; t++)
    {
        silargs[t].num_obs=num_obs;
        silargs[t].nmed=nmed;
        silargs[t].current_sil=&current_sil;
        silargs[t].nearest=&nearest;
        silargs[t].hist=&hist;
        silargs[t].silres=&silres;
        silargs[t].D = &D;
    }
    CreateAndRunThreadsWithDifferentArgs(nt,SilhoutteThread<disttype>,(void *)silargs,sizeof(SilhoutteThread_IO<disttype>));
    
    delete[] silargs;
 }
 Dt.EndClock(DEB & DEBPP); 

 Rcpp::NumericVector ret(num_obs);
 if (cl.hasAttribute("names"))
 {
  Rcpp::CharacterVector clnames=cl.names();
  if (clnames.length()>0)
  {
   if (clnames.length()!=ret.length())
    Rcpp::stop("Length of vector of names of the passed array of classifications is different from the length of the silhouette array.");
   ret.attr("names") = clnames;
  }
 }
 
 for (size_t t=0; t<current_sil.size(); t++)
  ret[t]=current_sil[t];

   
 return(ret);
}

//' CalculateSilhouette
//'
//' Calculates the silhouette of each point of those classified by a clustering algorithm.
//'
//' 
//' @param cl        The array of classification with the number of the class to which each point belongs to. This number must be in 1..number_of_classes.\cr
//'                  This function takes something like the L$clasif array which is the second element of the list returned by ApplyPAM
//' @param fdist     The binary file containing the symmetric matrix with the dissimilarities between cells (usually, generated by a call to CalcAndWriteDissimilarityMatrix)
//' @param nthreads  The number of used threads for parallel calculation.\cr
//'                  -1 means don't use threads (serial implementation).\cr
//'                   0 means let the program choose according to the number of cores and of points.\cr
//'                   Any other number forces this number of threads. Choosing more than the number of available cores is allowed, but discouraged.\cr
//'                   Default: 0
//' @return sil       Numeric vector with the values of the silhouette for each point, in the same order in which points are in cl.\cr
//'                   If cl is a named vector sil will be a named vector, too, with the same names.
//' @examples
//' # Synthetic problem: 10 random seeds with coordinates in [0..20]
//' # to which random values in [-0.1..0.1] are added
//' M<-matrix(0,100,500)
//' rownames(M)<-paste0("rn",c(1:100))
//' for (i in (1:10))
//' {
//'  p<-20*runif(500)
//'  Rf <- matrix(0.2*(runif(5000)-0.5),nrow=10)
//'  for (k in (1:10))
//'  {
//'   M[10*(i-1)+k,]=p+Rf[k,]
//'  }
//' }
//' tmpfile1=paste0(tempdir(),"/pamtest.bin")
//' JWriteBin(M,tmpfile1,dtype="float",dmtype="full")
//' tmpdisfile1=paste0(tempdir(),"/pamDL2.bin")
//' CalcAndWriteDissimilarityMatrix(tmpfile1,tmpdisfile1,distype="L2",restype="float",nthreads=0)
//' L <- ApplyPAM(tmpdisfile1,10,init_method="BUILD")
//' sil <- CalculateSilhouette(L$clasif,tmpdisfile1)
//' # Histogram of the silhouette. In this synthetic problem, almost 1 for all points
//' hist(sil)
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector CalculateSilhouette(Rcpp::NumericVector cl,std::string fdist,int nthreads=0)
{
 unsigned char mtype,ctype,e,md;
 indextype nr,nc;
 MatrixType(fdist,mtype,ctype,e,md,nr,nc);
 if (mtype!=MTYPESYMMETRIC)
  Rcpp::stop("This function can operate only with binary symmetric matrices.\n");
 if ((ctype!=FTYPE) && (ctype!=DTYPE))
  Rcpp::stop("This function can operate only with binary symmetric matrices with float or double elements.n");
  
 unsigned int nt=ChooseNumThreads(nthreads);
 
 if (ctype==FTYPE)
 {
  MemoryWarnings(nr,sizeof(float));
  return CalculateSilhouetteAux<float>(cl,fdist,nt);
 }
 else
 {
  MemoryWarnings(nr,sizeof(double));
  return CalculateSilhouetteAux<double>(cl,fdist,nt);
 }
}

//' NumSilToClusterSil
//'
//' Takes a silhouette in the form of a NumericVector, as returned by CalculateSilhouette, and returns it as a numeric matrix appropriate to be plotted by the package 'cluster'
//'
//' @param cl        The array of classification with the number of the class to which each point belongs to. This number must be in 1..number_of_classes.\cr
//'                  This function takes something like the L$clasif array which is the second element of the list returned by ApplyPAM
//' @param s         The numeric value of the silhouette for each point, with points in the same order as they appear in cl.\cr
//'                  This is the vector returned by a call to CalculateSilhouette with the same value of parameter cl.
//' @return sp       A silhouette in the format of the cluster package which is a NumericMatrix with as many rows as points and three columns: cluster, neighbor and sil_width.\cr
//'                  Its structure and dimension names are as in package 'cluster', which allows to use it with the silhouette plotting functions of such package\cr
//'                  This means you can do library(cluster) followed by plot(NumSilToClusterSil(cl,s)) to get a beatiful plot.
//' @examples
//' # Synthetic problem: 10 random seeds with coordinates in [0..20]
//' # to which random values in [-0.1..0.1] are added
//' M<-matrix(0,100,500)
//' rownames(M)<-paste0("rn",c(1:100))
//' for (i in (1:10))
//' {
//'  p<-20*runif(500)
//'  Rf <- matrix(0.2*(runif(5000)-0.5),nrow=10)
//'  for (k in (1:10))
//'  {
//'   M[10*(i-1)+k,]=p+Rf[k,]
//'  }
//' }
//' tmpfile1=paste0(tempdir(),"/pamtest.bin")
//' JWriteBin(M,tmpfile1,dtype="float",dmtype="full")
//' tmpdisfile1=paste0(tempdir(),"/pamDL2.bin")
//' CalcAndWriteDissimilarityMatrix(tmpfile1,tmpdisfile1,distype="L2",restype="float",nthreads=0)
//' L <- ApplyPAM(tmpdisfile1,10,init_method="BUILD")
//' sil <- CalculateSilhouette(L$clasif,tmpdisfile1)
//' sp <- NumSilToClusterSil(L$clasif,sil)
//' library(cluster)
//' plot(sp)
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix NumSilToClusterSil(Rcpp::NumericVector cl,Rcpp::NumericVector s)
{
 std::vector<silinfo> silres;
 
 for (indextype q=0; q<(unsigned int)s.length(); q++)
 {
  silinfo sl;
  sl.pnum=q+1;
  sl.ownclus=cl[q];
  sl.neiclus=0;
  sl.silvalue=s[q];
  silres.push_back(sl);
 }
 
 SortByClusterAndSilhouette(silres); 

 Rcpp::NumericMatrix retm(s.length(),3);

 Rcpp::CharacterVector cols(3);
 cols[0]="cluster";
 cols[1]="neighbor";
 cols[2]="sil_width";
 
 // The +1 is because the result must be returned in R-index convention
 for (indextype q=0;q<(unsigned int)s.length();q++)
 {
  retm(q,0)=silres[q].ownclus+1;
  retm(q,1)=silres[q].neiclus+1;
  retm(q,2)=silres[q].silvalue;
 }
 
 Rcpp::CharacterVector rows(s.length());
 for (indextype q=0;q<(unsigned int)s.length();q++)
  rows[q]=std::to_string(silres[q].pnum+1).c_str();
  
 Rcpp::List dimnames = Rcpp::List::create(rows,cols);
 retm.attr("dimnames")=dimnames;
 retm.attr("Ordered")=true;
 retm.attr("class")="silhouette";
 /* The next line was intended to create an object as similar as possible to the cluster::silhouette object
    of the cluster package but it provoked problems because the use of Rcpp::Language make internal
    calls to the R function SET_TYPEOF which is not in the R external API and this, in turn, provokes
    errors that prevent the package to be accepted by CRAN
 */
 //retm.attr("call")=Rcpp::Language("CalculateSilhouette","cl","fdist","nthreads");
 
 return(retm);
}

