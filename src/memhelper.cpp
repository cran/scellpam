#include<memhelper.h>

extern unsigned char DEB;

/*******************************************************
 * Function to get available memory
 * Uses R package memuse.
********************************************************/

// After call, avmem and avswap will return the available memory and available swap, as returned
// by the R package memuse through functions Sys.meminfo() and Sys.swapinfo(). Results are in KiB

void GetAvailableMemAndSwap(unsigned long &avmem,unsigned long &avswap)
{
 Rcpp::Function test("require");
 
 Rcpp::IntegerVector memuse_installed = test("memuse");
 if (DEB)                                                  // Messages on memory does not depend on the package purpose. They will be shown as long as any debug flag is active
 {
  Rcpp::Rcout << "Package memuse is ";
  if (memuse_installed[0]!=1)
   Rcpp::Rcout << "NOT installed. Cannot provide reliable memory information.\n";
  else
   Rcpp::Rcout << "installed. OK.\n";
 }
 if (memuse_installed[0]!=1)
 {
  avmem=avswap=0;
  Rcpp::warning("Package memuse if not installed. Cannot provide reliable memory information. Your request could exhaust your memory; not our fault. Install package 'memuse'.\n");
  return;
 } 
 
 Rcpp::Environment pkg = Rcpp::Environment::namespace_env("memuse");
 
 Rcpp::Function gm("Sys.meminfo");
 Rcpp::List gml=gm();
 Rcpp::S4 gml1=gml["freeram"];

 Rcpp::NumericVector frm=gml1.slot("size");
 std::string framunits=Rcpp::as<std::string>(gml1.slot("unit"));
 avmem=0;
 if (framunits=="GiB")
  avmem=frm[0]*1024*1024;
 if (framunits=="MiB")
  avmem=frm[0]*1024;
 if (framunits=="KiB")
  avmem=frm[0];
   
 Rcpp::Function gs("Sys.swapinfo");
 Rcpp::List gsl=gs();
 Rcpp::S4 gsl1=gsl["freeswap"];
 
 Rcpp::NumericVector fsw=gsl1.slot("size");
 std::string fswpunits=Rcpp::as<std::string>(gsl1.slot("unit"));
 avswap=0;
 if (fswpunits=="GiB")
  avswap=fsw[0]*1024*1024;
 if (fswpunits=="MiB")
  avswap=fsw[0]*1024;
 if (fswpunits=="KiB")
  avswap=fsw[0];
}

void MemoryWarnings(unsigned long nr,int s)
{
 unsigned long long estimated_size=s*nr*(nr+1)/2048;
 unsigned long mem=0,swap=0;
 
 GetAvailableMemAndSwap(mem,swap);
 
 if (mem==0)
  return;
  
 if (DEB)                        // Messages on memory does not depend on the package purpose. They will be shown as long as any debug flag is active
 {
  double percent=double(estimated_size)/double(mem);
  percent = int(10000.0*percent)/100.0;
  Rcpp::Rcout << "  Memory used by the matrix: " << estimated_size << " KiB, which is " << percent << "% of the available memory, which is " << mem << " Kib.\n";
  if (percent < 50.0)
   Rcpp::Rcout << "  That seems OK.\n";
  else
   if (percent < 75.0)
    Rcpp::Rcout << "  This is quite tight. Consider closing some application you don't need just now.\n";
   else
    Rcpp::Rcout << "  You are exhausting your memory. You should close some application you don't need just now.\n";
 }
 
 if (double(estimated_size)>double(mem)+double(swap))
  Rcpp::stop("Sorry, your computer has not enough memory to hold the matrix, not even using swap. Unfortunately, nothing can be done about that except getting more RAM.\n");
 
 if (double(estimated_size)>double(mem))
  Rcpp::warning("Your computer has not enough memory to hold the matrix so swap will be used. This means that calculation can be terribly slow. Use Ctrl-C to interrupt the program if you want.\n");
 
 if (double(estimated_size)>0.75*double(mem))
  Rcpp::warning("The matrix needs more than 75% of your computer's memory. This might provoke use of swap which will make calculation terribly slow. Close other applications, if possible, or interrupt the program with Ctrl-C.\n");
    
}

