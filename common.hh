#ifndef COMMON_H
#define COMMON_H

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
typedef blitz::Array<std::complex<double>,2> cmatrix;

typedef blitz::TinyVector<int,2> coord;
typedef blitz::Array<coord,1> coord_list;
typedef std::vector<coord> coord_vector;
typedef std::vector<coord_list> clist_vector;

template <typename T, int N_rank> inline
bool inrange(const blitz::Array<T,N_rank> &A, const blitz::TinyVector<int,N_rank> &I) 
{
	for (int i = 0; i < N_rank; i++)
		if ((I(i) < A.lbound(i)) || (I(i) > A.ubound(i))) return false;
	return true;
}

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

inline
unsigned int splitext(const std::string &in, std::string &froot, std::string &ext) {
	unsigned int fidx = in.rfind(".");
	if (fidx==std::string::npos)
		return -1;
	froot.assign(in.substr(0,fidx));
	ext.assign(in.substr(fidx));
	return fidx;
}
	

// template <typename T> inline
// T& MAX(T &x, T &y) {
// 	return (x > y) ? x : y;
// }

#endif
