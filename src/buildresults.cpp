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

#include <debugpar.h>
#include <vector>

extern unsigned char DEB;

/******************************************************************************
 * 
 * Functions related with build structures to facilitate statistical analysis
 *
 ***************************************************************************/

//' BuildAbundanceMatrix
//' 
//' Builds and returns a R matrix with as many rows as clusters and as many columns as groups in the set of cells (individuals).
//' The entry at row r, column c is the number if individuals of group c which the classifier has identified as belonging to cluster r
//'
//' @param clasif	The vector with the number of the cluster each cell belongs to. Usually obtained as L$clasif, being L the object returned by ApplyPAM.
//'                     It MUST be a vector of integers with as many components as cells and values in (1..number_of_clusters). Obviously, it can be a named
//'                     vector but the group names are not used.
//' @param gr		A numeric vector with the group (of those designed for the assay) to which each cell belongs to.
//'                     Normally obtained with GetSeuratGroups if the assay is in Seurat format. Otherwise, you will have to provide it yourself. It MUST be
//'                     vector of integers with as many components as cells and values in (1..number_of_groups). Obviously, it can be a named vector but
//'                     the cell names are not used.
//' @param expgroups    The expected number of groups. If it is left to its default value (which is 0) the number of groups is infered from parameter grname
//'                     as the maximum value in it. Otherwise, the passed value is used. This parameter is to prevent the extinction of some groups due to
//'                     previous expurge or filtering but whose trail we want to keep, even they are currently empty.
//' @return		M(numclusters,numgroups)  A R matrix as many rows as clusters and as many columns as groups
//' @examples
//' # Sorry, we can't provide examples here since they require the application to a real problem
//' # and therefore the load of the Seurat or splatter packages. Please, look at example in the
//' # vignette of this package.
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix BuildAbundanceMatrix(Rcpp::NumericVector clasif,Rcpp::IntegerVector gr,int expgroups=0)
{
 if (clasif.length()!=gr.length())
  Rcpp::stop("Lengths of vectors of clustering classification and group membership have not the same length (which must be the number of cells).\n");

 int gmax=gr[0];
 int gmin=gr[0];
 
 for (long int i=0; i<gr.length(); i++)
 {
  if (gr[i]>gmax)
   gmax=gr[i];
  if (gr[i]<gmin)
   gmin=gr[i];
 }
 if ((gmin!=1) || (gmax==gmin))
  Rcpp::stop("Vector of group membership minimal value is not 1, or maximal value is 1.\n");
 
 if (expgroups!=0)
 {
  if (gmax>expgroups)
   Rcpp::warning("More groups found in vector or groups than the expected number. We will keep the groups in the vector.\n");
  else
   gmax=expgroups;
 }
 
 int cmax=clasif[0];
 int cmin=clasif[0];
 for (long int i=0; i<clasif.length(); i++)
 {
  if (clasif[i]>cmax)
   cmax=clasif[i];
  if (clasif[i]<cmin)
   cmin=clasif[i];
 }
 if ((cmin!=1) || (cmax==cmin))
  Rcpp::stop("Vector of cluster membership minimal value is not 1, or maximal value is 1.\n");
  
 if (DEB & DEBSC)
  Rcpp::Rcout << clasif.length() << " cells distributed in " << cmax << " clusters and belonging to " << gmax << " groups.\n";
 
 Rcpp::NumericMatrix Ab(cmax,gmax);
 for (int i=0;i<int(cmax);i++)
  for (int j=0;j<gmax;j++)
   Ab(i,j)=0;
   
 for (int i=0;i<clasif.length();i++)
  Ab(clasif[i]-1,gr[i]-1)++;
  
 return(Ab);
}

