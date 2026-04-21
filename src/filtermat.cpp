#include <filtermat.h>

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

template void FilterDissim<float>(std::string ifname,std::string ofname,unsigned char mdinfo,std::vector<bool> keep,std::string addc);
template void FilterDissim<double>(std::string ifname,std::string ofname,unsigned char mdinfo,std::vector<bool> keep,std::string addc);
template void FilterDissim<long double>(std::string ifname,std::string ofname,unsigned char mdinfo,std::vector<bool> keep,std::string addc);

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

template void FilterCounts<unsigned short>(std::string ifname,bool issparse,unsigned char mdinfo,std::string ofname,std::vector<bool> keep,std::string addc);
template void FilterCounts<unsigned int>(std::string ifname,bool issparse,unsigned char mdinfo,std::string ofname,std::vector<bool> keep,std::string addc);
template void FilterCounts<unsigned long>(std::string ifname,bool issparse,unsigned char mdinfo,std::string ofname,std::vector<bool> keep,std::string addc);
template void FilterCounts<float>(std::string ifname,bool issparse,unsigned char mdinfo,std::string ofname,std::vector<bool> keep,std::string addc);
template void FilterCounts<double>(std::string ifname,bool issparse,unsigned char mdinfo,std::string ofname,std::vector<bool> keep,std::string addc);
template void FilterCounts<long double>(std::string ifname,bool issparse,unsigned char mdinfo,std::string ofname,std::vector<bool> keep,std::string addc);

