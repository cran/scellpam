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

#ifndef __CLOSECASES_H
#define __CLOSECASES_H
#include <dissimmat.h>
#include <cstdio>

struct args_to_close_thread
{
 indextype initial_row;
 indextype final_row;
 SymmetricMatrix<float> *M;
 indextype *uq;
 std::string *method;
 float *value;
 FullMatrix<indextype> *Cl;
};

#endif

