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

//' @importFrom memuse Sys.meminfo Sys.swapinfo

// This is the only place where this variable is declared. It is global for the full package and should only be changed by ScellpamSetDebug
unsigned char DEB=NODEBUG;

//' ScellpamSetDebug
//'
//' Sets debugging in scellpam package to ON (with TRUE) or OFF (with FALSE) for several parts of it.\cr
//' On package load the default status is OFF.\cr
//' Setting debugging of any part to ON shows a message. Setting to OFF does not show anything (since debugging is OFF...)
//'
//' @param deb       boolean, TRUE to generate debug messages for the scellpam (biological part) of this package and FALSE to turn them off. Default: true
//' @param debparpam boolean, TRUE to generate debug messages for the parallel PAM part inside this package and FALSE to turn them off. Default: false
//' @param debjmat   boolean, TRUE to generate debug messages for the jmatrix part inside this package and FALSE to turn them off. Default: false
//' @return          No return value, called for side effects (modification of internal variables)
//' @examples
//' ScellpamSetDebug(TRUE,debparpam=FALSE,debjmat=FALSE)
//' ScellpamSetDebug(TRUE,debparpam=TRUE,debjmat=FALSE)
//' ScellpamSetDebug(TRUE,debparpam=TRUE,debjmat=TRUE)
//' @export
// [[Rcpp::export]]
void ScellpamSetDebug(bool deb = true,bool debparpam = false,bool debjmat = false)
{
 if (deb)
 {
  DEB |= DEBSC;
  Rcpp::Rcout << "Debugging for scellpam (biological part) of the package set to ON.\n";
 }
 else
  DEB &= (~DEBSC);
  
 if (debparpam)
 {
  DEB |= DEBPP;
  Rcpp::Rcout << "Debugging for parallelpam inside scellpam package set to ON.\n";
 }
 else
  DEB &= (~DEBPP);
  
 if (debjmat)
 {
  DEB |= DEBJM;
  Rcpp::Rcout << "Debugging for jmatrix inside scellpam package set to ON.\n";
 }
 else
  DEB &= (~DEBJM);
}

//' ScellpamGetDebug
//'
//' Obtains the current state of the debugging parameter. To be used by R functions which want to adjust their messages according to the general debugging settings.
//'
//' @return     A list with the following keys:
//' \itemize{
//'      \item deb       - Boolean with the current debugging value of the scellpam (biological part) of this package
//'      \item debparpam - Boolean with the current debugging value of the parallel PAM part inside this package
//'      \item debjmat   - Boolean with the current debugging value of the jmatrix part inside this package
//' }
//' @examples
//' d<-ScellpamGetDebug()
//' @export
// [[Rcpp::export]]
Rcpp::List ScellpamGetDebug()
{
 bool d1,d2,d3;
 d1 = (DEB & DEBSC);
 d2 = (DEB & DEBPP);
 d3 = (DEB & DEBJM);
 Rcpp::List ret=Rcpp::List::create(Rcpp::Named("deb")=d1,Rcpp::Named("debparpam")=d2,Rcpp::Named("debjmat")=d3);
 return(ret);
}
