#ifndef COMMON_H
#define COMMON_H
/**
 * @file   common.hh
 * @author Daniel Meliza <dmeliza@uchicago.edu>
 * @date   Mon Mar  1 13:33:47 2010
 * 
 * @brief  Define common types and functions for znote
 * 
 * Copyright C Daniel Meliza, Z Chi 2010.  Licensed for use under Creative
 * Commons Attribution-Noncommercial-Share Alike 3.0 United States
 * License (http://creativecommons.org/licenses/by-nc-sa/3.0/us/).
 */
#include <blitz/array.h>
#include <cmath>
#include <complex>
#include <vector>

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327
#endif


typedef blitz::Array<double,2> dmatrix;
typedef blitz::Array<double,1> dvector;
typedef blitz::Array<int,1> ivector;
typedef blitz::Array<int,2> imatrix;
typedef blitz::Array<std::complex<double>,1> cvector;
typedef blitz::Array<std::complex<double>,2> cmatrix;

typedef blitz::TinyVector<int,2> coord;
typedef blitz::Array<coord,1> coord_list;
typedef std::vector<coord> coord_vector;
typedef std::vector<coord_list> clist_vector;
typedef std::vector<coord_vector> cvec_vector;

/** 
 * Determine if coordinates are in range for an array.
 * 
 * @param A  Array to check
 * @param I  Indices to check
 * 
 * @return true if A.lbound(i) <= I(i) <= A.ubound(i) for i in [0,N_rank)
 */
template <typename T, int N_rank> inline
bool inrange(const blitz::Array<T,N_rank> &A, const blitz::TinyVector<int,N_rank> &I) 
{
	for (int i = 0; i < N_rank; i++)
		if ((I(i) < A.lbound(i)) || (I(i) > A.ubound(i))) return false;
	return true;
}

/** 
 * Construct an evenly spaced vector across a range.
 * 
 * @param output Output vector. Resized and overwritten.
 * @param start Start point
 * @param stop Stop point (exclusive)
 * @param step Step size
 */
template <typename T> inline
void arange(blitz::Array<T,1> &output, int start, int stop, int step=1) 
{
	int n = (stop - start) / step;
	output.resize(n);
	output = start + blitz::tensor::i*step;
}

// template <class InputIterator>
// std::vector<std::iterator_traits<InputIterator>::value_type> 
// unique(InputIterator start, InputIterator stop, std::vector<std::iterator_traits<InputIterator>::value_type>)
// {
// 	std::iterator_traits<InputIterator>::value_type max_val = ;
// 	std::vector<std::iterator_traits<InputIterator>::value_type> out;
// 	for (++start;start!=stop;++start) {
// 		if (*start > max_val)

/** 
 * Split a file name into root and extension.
 * 
 * @param in input string
 * @param froot file name root
 * @param ext file name extension
 * 
 * @return split point
 */
inline
size_t splitext(const std::string &in, std::string &froot, std::string &ext) {
	size_t fidx = in.rfind(".");
	if (fidx==std::string::npos)
		return -1;
	froot.assign(in.substr(0,fidx));
	ext.assign(in.substr(fidx));
	return fidx;
}
	

#endif
