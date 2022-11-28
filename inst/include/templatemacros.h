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

/*
  This is an attempt to automatize to some extent the process of instantiation of tamplated constructors, destructors, operators and functions
  so that all the possible instances of the template are declared, and consequently compiled to object code to be put into a (shared) library.
  In that way any code will be able to use our library with any data type without recompilation or sorce code exposition.
  This is not because we don't want our code to be exposed; it is open-source, anyway. The problem is that Rcpp compiles in a faster and safest
  way. The inconvenient is the size of the generated .so library but it is not so big, and the savings in compile time is probably worth.
  */
#ifndef _TEMPLATE_MACROS_H
#define _TEMPLATE_MACROS_H

#define SINGLE_ARG(...) __VA_ARGS__

#define TEMPLATES_CONST(B,C) \
template B<unsigned char>::B(C); \
template B<char>::B(C); \
template B<unsigned short>::B(C); \
template B<short>::B(C); \
template B<unsigned int>::B(C); \
template B<int>::B(C); \
template B<unsigned long>::B(C); \
template B<long>::B(C); \
template B<float>::B(C); \
template B<double>::B(C); \
template B<long double>::B(C);

#define TEMPLATES_DEFAULT_DEST(B) \
template B<unsigned char>::~B(); \
template B<char>::~B(); \
template B<unsigned short>::~B(); \
template B<short>::~B(); \
template B<unsigned int>::~B(); \
template B<int>::~B(); \
template B<unsigned long>::~B(); \
template B<long>::~B(); \
template B<float>::~B(); \
template B<double>::~B(); \
template B<long double>::~B();

#define TEMPLATES_COPY_CONST(A) \
template A<unsigned char>::A( const A<unsigned char>& other); \
template A<char>::A( const A<char>& other); \
template A<unsigned short>::A( const A<unsigned short>& other); \
template A<short>::A( const A<short>& other); \
template A<unsigned int>::A( const A<unsigned int>& other); \
template A<int>::A( const A<int>& other); \
template A<unsigned long>::A( const A<unsigned long>& other); \
template A<long>::A( const A<long>& other); \
template A<float>::A( const A<float>& other); \
template A<double>::A( const A<double>& other); \
template A<long double>::A( const A<long double>& other);

#define TEMPLATES_OPERATOR(A,B) \
template A<unsigned char>& A<unsigned char>::operator B( const A<unsigned char>& other); \
template A<char>& A<char>::operator B( const A<char>& other); \
template A<unsigned short>& A<unsigned short>::operator B( const A<unsigned short>& other); \
template A<short>& A<short>::operator B( const A<short>& other); \
template A<unsigned int>& A<unsigned int>::operator B( const A<unsigned int>& other); \
template A<int>& A<int>::operator B( const A<int>& other); \
template A<unsigned long>& A<unsigned long>::operator B( const A<unsigned long>& other); \
template A<long>& A<long>::operator B( const A<long>& other); \
template A<float>& A<float>::operator B( const A<float>& other); \
template A<double>& A<double>::operator B( const A<double>& other); \
template A<long double>& A<long double>::operator B( const A<long double>& other);

#define TEMPLATES_FUNC(A,B,C,D) \
template A B<unsigned char>::C(D); \
template A B<char>::C(D); \
template A B<unsigned short>::C(D); \
template A B<short>::C(D); \
template A B<unsigned int>::C(D); \
template A B<int>::C(D); \
template A B<unsigned long>::C(D); \
template A B<long>::C(D); \
template A B<float>::C(D); \
template A B<double>::C(D); \
template A B<long double>::C(D);

#define TEMPLATES_SETFUNC(A,B,C,D,E) \
template A B<unsigned char>::C(D,unsigned char E); \
template A B<char>::C(D,char E); \
template A B<unsigned short>::C(D,unsigned short E); \
template A B<short>::C(D,short E); \
template A B<unsigned int>::C(D,unsigned int E); \
template A B<int>::C(D,int E); \
template A B<unsigned long>::C(D,unsigned long E); \
template A B<long>::C(D,long E); \
template A B<float>::C(D,float E); \
template A B<double>::C(D,double E); \
template A B<long double>::C(D,long double E);

#define TEMPLATES_SETFUNCVEC(A,B,C,D,E) \
template A B<unsigned char>::C(D,std::vector<unsigned char> E); \
template A B<char>::C(D,std::vector<char> E); \
template A B<unsigned short>::C(D,std::vector<unsigned short> E); \
template A B<short>::C(D,std::vector<short> E); \
template A B<unsigned int>::C(D,std::vector<unsigned int> E); \
template A B<int>::C(D,std::vector<int> E); \
template A B<unsigned long>::C(D,std::vector<unsigned long> E); \
template A B<long>::C(D,std::vector<long> E); \
template A B<float>::C(D,std::vector<float> E); \
template A B<double>::C(D,std::vector<double> E); \
template A B<long double>::C(D,std::vector<long double> E);

#define TEMPLATES_FUNCR(B,C,D) \
template unsigned char B<unsigned char>::C(D); \
template char B<char>::C(D); \
template unsigned short B<unsigned short>::C(D); \
template short B<short>::C(D); \
template unsigned int B<unsigned int>::C(D); \
template int B<int>::C(D); \
template unsigned long B<unsigned long>::C(D); \
template long B<long>::C(D); \
template float B<float>::C(D); \
template double B<double>::C(D); \
template long double B<long double>::C(D);

#define TEMPLATES_FUNCRCONST(B,C,D) \
template unsigned char B<unsigned char>::C(D) const; \
template char B<char>::C(D) const; \
template unsigned short B<unsigned short>::C(D) const; \
template short B<short>::C(D) const; \
template unsigned int B<unsigned int>::C(D) const; \
template int B<int>::C(D) const; \
template unsigned long B<unsigned long>::C(D) const; \
template long B<long>::C(D) const; \
template float B<float>::C(D) const; \
template double B<double>::C(D) const; \
template long double B<long double>::C(D) const;
#endif

