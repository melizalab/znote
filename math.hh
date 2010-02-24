#ifndef MATH_H
#define MATH_H

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
typedef blitz::TinyVector<int,2> coord;
typedef blitz::Array<coord,1> coord_list;
typedef std::vector<coord> coord_vector;
typedef std::vector<coord_list> clist_vector;

template <typename T, int N_rank> inline
bool inrange(const blitz::Array<T,N_rank> &A, const blitz::TinyVector<int,N_rank> &I) {
	for (int i = 0; i < N_rank; i++)
		if ((I(i) < A.lbound(i)) || (I(i) > A.ubound(i))) return false;
	return true;
}

template <typename T> inline
T& MAX(T &x, T &y) {
	return (x > y) ? x : y;
}

// template <typename T1, typename T2>
// bool same(const T1 &v1, const T2 &v2) {
// 	return v1==v2;
// }

// template <typename T1, typename T2, int N_rank> inline
// bool same(const blitz::TinyVector<T1,N_rank> &A1, const blitz::TinyVector<T2,N_rank> &A2) {
// 	for (int i = 0; i < N_rank; i++)
// 		if (A1(i)!=A2(i)) return false;
// 	return true;
// }

#endif
