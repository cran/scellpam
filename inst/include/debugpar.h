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

#ifndef _DEBUGPAR_H
#define _DEBUGPAR_H

#include  <Rcpp.h>

// These are constants to allow selective debug by sub-package.
// Each package will print messages or not using a test with logical AND
// between its particular constant and the DEB global variable. This
// allows the use of the system either in each separate package or
// in the global one
const unsigned char NODEBUG=0x0;
const unsigned char DEBJM=0x01;
const unsigned char DEBPP=0x02;
const unsigned char DEBSC=0x04;

#endif
