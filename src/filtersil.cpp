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

extern unsigned char DEB;

// This function returns the R object with the modified L$med and L$clasif fields and updates a vector with the _indexes_ of the points
// which must be kept, not with its values (which we don't care about)
Rcpp::List FilterByQuantile(Rcpp::NumericVector s, float q, Rcpp::List L,std::vector<bool> &kept)
{
 float threshold;
 {
  Rcpp::NumericVector iq = Rcpp::clone(s);
  iq.sort();
  threshold=iq(indextype(0.5+q*float(s.length())));
 }
 
 indextype remain=0;
 for (indextype c=0; c<indextype(s.length()); c++)
  if (s(c)>=threshold)
  {
   kept[c]=true;
   remain++;
  }
  
 // No need to clone here, we do not change L$med. We will return a new vector
 Rcpp::NumericVector med = L["med"];
 Rcpp::StringVector mnames = med.names();
 
 Rcpp::NumericVector clasif = L["clasif"];
 Rcpp::StringVector cnames = clasif.names();
 
 if (clasif.length() != s.length())
  Rcpp::stop("Length of silhouette and of vector of classified points in field $clasif are not equal. Are you sure silhouette corresponds to this clustering?\n");
  
 // This is a filter to be sure we don't get rid of any medoid.
 indextype nfmed=0;
 for (indextype t=0; t<indextype(med.length()); t++)
  if (kept[indextype(med(t)-1)]==false)
  {
   nfmed++;
   kept[indextype(med(t)-1)]=true;
  }
 
 remain=0;
 for (indextype c=0; c<indextype(s.length()); c++)
  if (kept[c])
   remain++;
   
 if (DEB & DEBSC)
  Rcpp::Rcout << "After filtering silhouette with quantile " << q << " (threshold " << threshold << ") " << remain << " of the " << s.length() << " points remain.\n";
 
 if (nfmed!=0)
 {
  std::ostringstream errst;
  errst << nfmed << " of the medoids have been kept, even they were below the threshold (which seems problematic. Check your clusters...).\n";
  Rcpp::warning(errst.str());
  if (DEB & DEBSC)
   Rcpp::Rcout << nfmed << " of the medoids have been kept, even they were below the threshold (which seems problematic. Check your clusters...).\n";
 }
 
 // The clasif is simply copied, getting rid of the removed points.
 // The points do not change their assignment to a cluster, and clasif contains cluster numbers (identifiers) in R-index, as they must be in the returned variable.
 
 Rcpp::NumericVector new_clasif(remain);
 // The names of the cells must be known, to keep those of the cells that must be kept.
 Rcpp::StringVector new_clasif_names(remain);
 
 Rcpp::NumericVector new_medoid_indices(med.length()); 
 indextype push_index=0;
 indextype m1;
 for (indextype t=0; t<kept.size(); t++)
 {
  if (kept[t])
  {
   if (push_index>remain)
    Rcpp::stop("Too many points kept..??? (unexpected error)\n");
   
   new_clasif(push_index)=clasif(t);  // We do not add/substract 1 here, since in both vectors indices are R-type (from 1)
   if (cnames.length()>0)
    new_clasif_names(push_index)=cnames(t);
    
   // Where in the new clasif have medoids been left?
   m1=0;
   while (m1<indextype(med.length()) && med(m1)-1!=t)
    m1++;
   if (m1<indextype(med.length()))   // Point t is a medoid, it has been found...
    new_medoid_indices(m1)=push_index+1;
   push_index++;
  }
  
 }
  
 if (cnames.length()>0) 
  new_clasif.names() = new_clasif_names;
 
 if (mnames.length()>0)
  new_medoid_indices.names() = mnames;
  
 Rcpp::List ret;
 
 ret["med"]=new_medoid_indices;
 ret["clasif"] = new_clasif;
 
 return ret;
}

// This function returns the R object with the modified L$med and L$clasif fields and updates a vector with the _indexes_ of the points
// which must be kept, not with its values (which we don't care about)
Rcpp::List FilterByThreshold(Rcpp::NumericVector s, float threshold, Rcpp::List L,std::vector<bool> &kept)
{
 indextype remain=0;
 for (indextype c=0; c<indextype(s.length()); c++)
  if (s(c)>=threshold)
  {
   kept[c]=true;
   remain++;
  }
  
 // No need to clone here, we do not change L$med. We will return a new vector
 Rcpp::NumericVector med = L["med"];
 Rcpp::StringVector mnames = med.names();
 
 Rcpp::NumericVector clasif = L["clasif"];
 Rcpp::StringVector cnames = clasif.names();
 
 if (clasif.length() != s.length())
  Rcpp::stop("Length of silhouette and of vector of classified points in field $clasif are not equal. Are you sure silhouette corresponds to this clustering?\n");
  
 // This is a filter to be sure we don't get rid of any medoid.
 indextype nfmed=0;
 for (indextype t=0; t<indextype(med.length()); t++)
  if (kept[indextype(med(t)-1)]==false)
  {
   nfmed++;
   kept[indextype(med(t)-1)]=true;
  }
 
 remain=0;
 for (indextype c=0; c<indextype(s.length()); c++)
  if (kept[c])
   remain++;
   
 if (DEB & DEBSC)
  Rcpp::Rcout << "After filtering silhouette with threshold " << threshold << ", " << remain << " of the " << s.length() << " points remain.\n";
 
 if (nfmed!=0)
 {
  std::ostringstream errst;
  errst << nfmed << " of the medoids have been kept, even they were below the threshold (which seems problematic. Check your clusters...).\n";
  Rcpp::warning(errst.str());
  if (DEB & DEBSC)
   Rcpp::Rcout << nfmed << " of the medoids have been kept, even they were below the threshold (which seems problematic. Check your clusters...).\n";
 }
 
 // The clasif is simply copied, getting rid of the removed points.
 // The points do not change their assignment to a cluster, and clasif contains cluster numbers (identifiers) in R-index, as they must be in the returned variable.
 
 Rcpp::NumericVector new_clasif(remain);
 // The names of the cells must be known, to keep those of the cells that must be kept.
 Rcpp::StringVector new_clasif_names(remain);
 
 Rcpp::NumericVector new_medoid_indices(med.length()); 
 indextype push_index=0;
 indextype m1;
 for (indextype t=0; t<kept.size(); t++)
 {
  if (kept[t])
  {
   if (push_index>remain)
    Rcpp::stop("Too many points kept..??? (unexpected error)\n");
   
   new_clasif(push_index)=clasif(t);  // We do not add/substract 1 here, since in both vectors indices are R-type (from 1)
   if (cnames.length()>0)
    new_clasif_names(push_index)=cnames(t);
    
   // Where in the new clasif have medoids been left?
   m1=0;
   while (m1<indextype(med.length()) && med(m1)-1!=t)
    m1++;
   if (m1<indextype(med.length()))   // Point t is a medoid, it has been found...
    new_medoid_indices(m1)=push_index+1;
   push_index++;
  }
  
 }
  
 if (cnames.length()>0) 
  new_clasif.names() = new_clasif_names;
 
 if (mnames.length()>0)
  new_medoid_indices.names() = mnames;
  
 Rcpp::List ret;
 
 ret["med"]=new_medoid_indices;
 ret["clasif"] = new_clasif;
 
 return ret;
}

template <typename T>
void FilterCounts(std::string ifname,bool issparse,unsigned char mdinfo,std::string ofname,std::vector<bool> keep,std::string addc)
{
 indextype nrt,nrf,nc,nf;
 
 nrf=0;
 for (indextype i=0;i<keep.size();i++)
  if (keep[i])
   nrf++;
   
 if (issparse)
 {
  SparseMatrix<T> M(ifname);
  nrt=M.GetNRows();
  nc=M.GetNCols();
  
  if (nrt==nrf)
  {
   // WARNING: Check here. Somthing in the copy constructor???
   //SparseMatrix<T> Mf(M);
   //Mf.WriteBin(ofname);
   if (addc != "")
   {
    if (mdinfo & COMMENT)
     M.SetComment(M.GetComment()+addc);
    else
     M.SetComment(addc);
   }
   M.WriteBin(ofname);
  }
  else
  { 
   SparseMatrix<T> Mf(nrf,nc);
   nf=0;
   for (indextype r=0; r<nrt; r++)
    if (keep[r])
    {
     for (indextype c=0; c<nc; c++)
      Mf.Set(nf,c,M.Get(r,c));
     nf++;
    }
   
   if (mdinfo & COL_NAMES)
    Mf.SetColNames(M.GetColNames());
   if (mdinfo & ROW_NAMES)
   {
    std::vector<std::string> names=M.GetRowNames();
    std::vector<std::string> rem_cells;
    for (indextype r=0; r<nrt; r++)
     if (keep[r])
      rem_cells.push_back(names[r]);
    Mf.SetRowNames(rem_cells);
   }
   if (mdinfo & COMMENT)
    Mf.SetComment(M.GetComment()+addc);
   else
    if (addc!="")
     Mf.SetComment(addc);
     
   Mf.WriteBin(ofname);
  }
 }
 else
 {
  FullMatrix<T> M(ifname);
  nrt=M.GetNRows();
  nc=M.GetNCols();
  
  if (nrt==nrf)
  {
   // WARNING: Check here. Somthing in the copy constructor???
   //FullMatrix<T> Mf(M);
   //Mf.WriteBin(ofname);
   if (addc != "")
   {
    if (mdinfo & COMMENT)
     M.SetComment(M.GetComment()+addc);
    else
     M.SetComment(addc);
   }
   M.WriteBin(ofname);
  }
  else
  {
   FullMatrix<T> Mf(nrf,nc);
   nf=0;
   for (indextype r=0; r<nrt; r++)
    if (keep[r])
    {
     for (indextype c=0; c<nc; c++)
      Mf.Set(nf,c,M.Get(r,c));
     nf++;
    }
   if (mdinfo & COL_NAMES)
    Mf.SetColNames(M.GetColNames());
   if (mdinfo & ROW_NAMES)
   {
    std::vector<std::string> names=M.GetRowNames();
    std::vector<std::string> rem_cells;
    for (indextype r=0; r<nrt; r++)
     if (keep[r])
      rem_cells.push_back(names[r]);
    Mf.SetRowNames(rem_cells);
   }
   if (mdinfo & COMMENT)
    Mf.SetComment(M.GetComment()+addc);
   else
    if (addc!="")
     Mf.SetComment(addc);
      
   Mf.WriteBin(ofname);
  }
 }
}

template <typename T>
void FilterDissim(std::string ifname,std::string ofname,unsigned char mdinfo,std::vector<bool> keep,std::string addc)
{
 SymmetricMatrix<T> M(ifname);
 indextype nrt=M.GetNRows();
 
 indextype nrf=0;
 for (indextype i=0; i<keep.size(); i++)
  if (keep[i])
   nrf++;
   
 SymmetricMatrix<T> Mf(nrf);
 
 indextype nf,nc;
 nf=0;
 for (indextype r=0; r<nrt; r++)
  if (keep[r])
  {
   nc=0;
   for (indextype c=0; c<=r; c++)
    if (keep[c])
    {
     Mf.Set(nf,nc,M.Get(r,c));
     nc++;
    }
   nf++;
  } 
   
 if (mdinfo & ROW_NAMES)
 {
  std::vector<std::string> names=M.GetRowNames();
  std::vector<std::string> rem_cells;  
  for (indextype r=0; r<nrt; r++)
   if (keep[r])
    rem_cells.push_back(names[r]);
  Mf.SetRowNames(rem_cells);
 }
 if (mdinfo & COMMENT)
  Mf.SetComment(M.GetComment()+addc); 
 else
  if (addc!="")
   Mf.SetComment(addc);
  
 Mf.WriteBin(ofname);  
}

//' FilterBySilhouetteQuantile
//'
//' Takes a silhouette, as returned by CalculateSilhouette, the list of medoids and class assignments, as returned by ApplyPam,
//' a quantile and the matrices of counts and dissimilarities and constructs the corresponding matrices clearing off the points (cells) whose silhoutte is
//' below the lower quantile, except if they are medoids.\cr
//'
//' The renumbering of indices in the returned cluster may seem confusing at first but it was the way of fitting this with the rest
//' of the package. Anyway, notice that if the numeric vectors in the input parameter L were named vectors, the cells names are appropriately kept
//' in the result so cell identity is preserved. Moreover, if the counts and dissimilarity input matrices had row and/or column names, they
//' are preserved in the filtered matrices, too.
//' 
//' @param s          A numeric vector with the sihouette coefficient of each point (cell) in a classification, as returned by CalculateSilhouette.
//' @param L          A list of two numeric vectors, L$med and L$clasif, obtained normally as the object returned by ApplyPAM.
//' @param fallcounts A string with the name of the binary file containing the matrix of counts per cell. It can be either a full or a sparse matrix.
//' @param ffilcounts A string with the name of the binary file that will contain the selected cells. It will have the same character (full/sparse) and type of the complete file.
//' @param falldissim A string with the name of the binary file containing the dissimilarity matrix of the complete set of cells. It must be a symmetric matrix of floats.
//' @param ffildissim A string with the name of the binary file that will contain  the dissimilarity matrix for the remaining cells. It will be a symmetric matrix of floats.
//' @param q          Quantile to filter. All points (cells) whose silhouette is below this quantile will be filtered out. Default: 0.2
//' @param addcom     Boolean to indicate if a comment must be appended to the current comment of counts and dissimilarity matrices to indicate that they are the result of a filtering process. This comment is automatically generated and contains the value of quantile q. Succesive applications add comments at the end of those already present. Default: TRUE
//' @return Lr["med","clasif"] A list of two numeric vectors.\cr
//'                      Lr$med is a modification of the correponding first element of the passed L parameter.\cr
//'                      Lr$clasif has as many components as remaining instances.\cr
//'                      Since points (cells) will have been removed, medoid numbering is modified. Therefore, Lr$med has the NEW index of each medoid in the filtered set.\cr
//'                      Lr$clasif contains the number of the medoid (i.e.: the cluster) to which each instance has been assigned, and therefore does not change.\cr
//'                      All indexes start at 1 (R convention). Please, see Details section\cr
//'
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
//' tmpdisfile1=paste0(tempdir(),"/pamDl2.bin")
//' CalcAndWriteDissimilarityMatrix(tmpfile1,tmpdisfile1,distype="L2",restype="float",nthreads=0)
//' L <- ApplyPAM(tmpdisfile1,10,init_method="BUILD")
//' # Which are the medoids
//' L$med
//' sil <- CalculateSilhouette(L$clasif,tmpdisfile1)
//' tmpfiltfile1=paste0(tempdir(),"/pamtestfilt.bin")
//' tmpfiltdisfile1=paste0(tempdir(),"/pamDL2filt.bin")
//' Lf<-FilterBySilhouetteQuantile(sil,L,tmpfile1,tmpfiltfile1,tmpdisfile1,tmpfiltdisfile1,
//'                                q=0.4,addcom=TRUE)
//' # The new medoids are the same points but renumbered, since the L$clasif array has less points
//' Lf$med
//' @export
// [[Rcpp::export]]
Rcpp::List FilterBySilhouetteQuantile(Rcpp::NumericVector s, Rcpp::List L,std::string fallcounts, std::string ffilcounts, std::string falldissim, std::string ffildissim, float q=0.2,bool addcom=true)
{   
 if (q<=0.0 || q>=1.0)
  Rcpp::stop("quantile must be a floating-point number 0.0 < q < 1.0.\n");
  
 std::vector<bool> selected(s.length(),false);
 
 Rcpp::List ret=FilterByQuantile(s,q,L,selected);

 // Characteristics of the counts matrix
 unsigned char mtc,ctc,endianc,mdinfoc;
 indextype nrc,ncc;
 MatrixType(fallcounts,mtc,ctc,endianc,mdinfoc,nrc,ncc);
 
 // Characteristics of the dissimilarity matrix
 unsigned char mtd,ctd,endiand,mdinfod;
 indextype nrd,ncd;
 MatrixType(falldissim,mtd,ctd,endiand,mdinfod,nrd,ncd);
 
 if (mtc!=MTYPESPARSE && mtc!=MTYPEFULL)
  Rcpp::stop("The count matrix contained in the input counts file is neither full nor sparse. Not a valid count matrix.\n");
  
 if (mtd!=MTYPESYMMETRIC)
  Rcpp::stop("The count matrix contained in the input dissimilarity file is not a symmetric matrix. Not a valid distance/dissimilarity matrix.\n");
  
 if (nrc!=nrd)
  Rcpp::stop("The number of points in the count matrix and that in the dissimilarity matrix are different. Are you sure they corespond to the same data?.\n");
  
 if (nrc!=(unsigned int)s.length())
  Rcpp::stop("The number of points in the silhouette vector and the number of rows in the count matrix are not equal.\nWas this silhouette calculated with such points?.\n");  
 
 std::string added_comment="";
 if (addcom)
 {
  std::ostringstream stcom;
  stcom << " Filtered by silhouette from file " << fallcounts << " with quantile " << q << ". ";
  added_comment = stcom.str();
 }
 switch (ctc)
 {
  case USTYPE: FilterCounts<unsigned short>(fallcounts,(mtc==MTYPESPARSE),mdinfoc,ffilcounts,selected,added_comment); break;
  case UITYPE: FilterCounts<unsigned int>(fallcounts,(mtc==MTYPESPARSE),mdinfoc,ffilcounts,selected,added_comment); break;
  case ULTYPE: FilterCounts<unsigned long>(fallcounts,(mtc==MTYPESPARSE),mdinfoc,ffilcounts,selected,added_comment); break;
  case FTYPE:  FilterCounts<float>(fallcounts,(mtc==MTYPESPARSE),mdinfoc,ffilcounts,selected,added_comment); break;
  case DTYPE:  FilterCounts<double>(fallcounts,(mtc==MTYPESPARSE),mdinfoc,ffilcounts,selected,added_comment); break;
  case LDTYPE:  FilterCounts<long double>(fallcounts,(mtc==MTYPESPARSE),mdinfoc,ffilcounts,selected,added_comment); break;
  default: Rcpp::stop("The count matrix is not of any of the valid data types (unsigned short, unsigned int, unsigned long, float, double or longdouble).\n"); break;
 }


 if (addcom)
 {
  std::ostringstream stcom;
  stcom << " Filtered by silhouette from file " << falldissim << " with quantile " << q << ". ";
  added_comment = stcom.str();
 }
 switch (ctd)
 {
  case FTYPE: FilterDissim<float>(falldissim,ffildissim,mdinfod,selected,added_comment); break;
  case DTYPE: FilterDissim<double>(falldissim,ffildissim,mdinfod,selected,added_comment); break; 
  case LDTYPE: FilterDissim<long double>(falldissim,ffildissim,mdinfod,selected,added_comment); break; 
  default: Rcpp::stop("The dissimilarity matrix is not of any of the valid data types (float, double or longdouble).\n"); break;
 }
 
 return(ret);
}

//' FilterBySilhouetteThreshold
//'
//' Takes a silhouette, as returned by CalculateSilhouette, the list of medoids and class assignments, as returned by ApplyPam,
//' a threshold and the matrices of counts and dissimilarities and constructs the corresponding matrices clearing off the points (cells) whose silhoutte is
//' below the threshold, except if they are medoids.\cr
//'
//' The renumbering of indices in the returned cluster may seem confusing at first but it was the way of fitting this with the rest
//' of the package. Anyway, notice that if the numeric vectors in the input parameter L were named vectors, the cells names are appropriately kept
//' in the result so cell identity is preserved. Moreover, if the counts and dissimilarity input matrices had row and/or column names, they
//' are preserved in the filtered matrices, too.
//' 
//' @param s          A numeric vector with the sihouette coefficient of each point in a classification, as returned by CalculateSilhouette.
//' @param L          A list of two numeric vectors, L$med and L$clasif, obtained normally as the object returned by ApplyPAM.
//' @param fallcounts A string with the name of the binary file containing the matrix of counts per cell. It can be either a full or a sparse matrix.
//' @param ffilcounts A string with the name of the binary file that will contain the selected cells. It will have the same character (full/sparse) and type of the complete file.
//' @param falldissim A string with the name of the binary file containing the dissimilarity matrix of the complete set of cells. It must be a symmetric matrix of floats.
//' @param ffildissim A string with the name of the binary file that will contain  the dissimilarity matrix for the remaining cells. It will be a symmetric matrix of floats.
//' @param thres      Threshold to filter. All points whose silhouette is below this threshold will be filtered out. Default: 0.0 (remember that silhouette is in [-1..1])
//' @param addcom     Boolean to indicate if a comment must be appended to the current comment of counts and dissimilarity matrices to indicate that they are the result of a filtering process. This comment is automatically generated and contains the value of threshold t. Succesive applications add comments at the end of those already present. Default: TRUE
//' @return Lr["med","clasif"] A list of two numeric vectors.\cr
//'                      Lr$med is a modification of the correponding first element of the passed L parameter.\cr
//'                      Lr$clasif has as many components as remaining instances.\cr
//'                      Since points will have been removed, medoid numbering is modified. Therefore, Lr$med has the NEW index of each medoid in the filtered set.\cr
//'                      Lr$clasif contains the number of the medoid (i.e.: the cluster) to which each instance has been assigned, and therefore does not change.\cr
//'                      All indexes start at 1 (R convention). Please, see Details section\cr
//'
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
//' tmpdisfile1=paste0(tempdir(),"/pamDl2.bin")
//' CalcAndWriteDissimilarityMatrix(tmpfile1,tmpdisfile1,distype="L2",restype="float",nthreads=0)
//' L <- ApplyPAM(tmpdisfile1,10,init_method="BUILD")
//' # Which are the medoids
//' L$med
//' sil <- CalculateSilhouette(L$clasif,tmpdisfile1)
//' tmpfiltfile1=paste0(tempdir(),"/pamtestfilt.bin")
//' tmpfiltdisfile1=paste0(tempdir(),"/pamDL2filt.bin")
//' Lf<-FilterBySilhouetteThreshold(sil,L,tmpfile1,tmpfiltfile1,tmpdisfile1,tmpfiltdisfile1,
//'                                thres=0.4,addcom=TRUE)
//' # The new medoids are the same points but renumbered, since the L$clasif array has less points
//' Lf$med
//' @export
// [[Rcpp::export]]
Rcpp::List FilterBySilhouetteThreshold(Rcpp::NumericVector s, Rcpp::List L,std::string fallcounts, std::string ffilcounts, std::string falldissim, std::string ffildissim, float thres=0.0,bool addcom=true)
{   
 if (thres<-1.0 || thres>1.0)
  Rcpp::stop("threshold must be a floating-point number -1.0 <= q <= 1.0.\n");
  
 std::vector<bool> selected(s.length(),false);
 
 Rcpp::List ret=FilterByThreshold(s,thres,L,selected);

 // Characteristics of the counts matrix
 unsigned char mtc,ctc,endianc,mdinfoc;
 indextype nrc,ncc;
 MatrixType(fallcounts,mtc,ctc,endianc,mdinfoc,nrc,ncc);
 
 // Characteristics of the dissimilarity matrix
 unsigned char mtd,ctd,endiand,mdinfod;
 indextype nrd,ncd;
 MatrixType(falldissim,mtd,ctd,endiand,mdinfod,nrd,ncd);
 
 if (mtc!=MTYPESPARSE && mtc!=MTYPEFULL)
  Rcpp::stop("The count matrix contained in the input counts file is neither full nor sparse. Not a valid count matrix.\n");
  
 if (mtd!=MTYPESYMMETRIC)
  Rcpp::stop("The count matrix contained in the input dissimilarity file is not a symmetric matrix. Not a valid distance/dissimilarity matrix.\n");
  
 if (nrc!=nrd)
  Rcpp::stop("The number of points in the count matrix and that in the dissimilarity matrix are different. Are you sure they corespond to the same data?.\n");
  
 if (nrc!=(unsigned int)s.length())
  Rcpp::stop("The number of points in the silhouette vector and the number of rows in the count matrix are not equal.\nWas this silhouette calculated with such points?.\n");  
 
 std::string added_comment="";
 if (addcom)
 {
  std::ostringstream stcom;
  stcom << " Filtered by silhouette from file " << fallcounts << " with threshold " << thres << ". ";
  added_comment = stcom.str();
 }
 
 switch (ctc)
 {
  case USTYPE: FilterCounts<unsigned short>(fallcounts,(mtc==MTYPESPARSE),mdinfoc,ffilcounts,selected,added_comment); break;
  case UITYPE: FilterCounts<unsigned int>(fallcounts,(mtc==MTYPESPARSE),mdinfoc,ffilcounts,selected,added_comment); break;
  case ULTYPE: FilterCounts<unsigned long>(fallcounts,(mtc==MTYPESPARSE),mdinfoc,ffilcounts,selected,added_comment); break;
  case FTYPE:  FilterCounts<float>(fallcounts,(mtc==MTYPESPARSE),mdinfoc,ffilcounts,selected,added_comment); break;
  case DTYPE:  FilterCounts<double>(fallcounts,(mtc==MTYPESPARSE),mdinfoc,ffilcounts,selected,added_comment); break;
  case LDTYPE:  FilterCounts<long double>(fallcounts,(mtc==MTYPESPARSE),mdinfoc,ffilcounts,selected,added_comment); break;
  default: Rcpp::stop("The count matrix is not of any of the valid data types (unsigned short, unsigned int, unsigned long, float, double or longdouble).\n"); break;
 }

 if (addcom)
 {
  std::ostringstream stcom;
  stcom << " Filtered by silhouette from file " << falldissim << " with threshold " << thres << ". ";
  added_comment = stcom.str();
 }
 
 switch (ctd)
 {
  case FTYPE: FilterDissim<float>(falldissim,ffildissim,mdinfod,selected,added_comment); break;
  case DTYPE: FilterDissim<double>(falldissim,ffildissim,mdinfod,selected,added_comment); break; 
  case LDTYPE: FilterDissim<long double>(falldissim,ffildissim,mdinfod,selected,added_comment); break; 
  default: Rcpp::stop("The dissimilarity matrix is not of any of the valid data types (float, double or longdouble).\n"); break;
 }
 
 return(ret);
}

