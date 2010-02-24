#include <cstdlib>
#include <cstdio>
#include <string.h>
#include <cmath>
#include <float.h>
#include <complex>
#include "spect.hh"

/* some LAPACK prototypes */
extern void dsterf_(int *N, double *D, double *E, int *INFO);

extern void dgtsv_(int *N, int *NRHS,
		   double *DL, double *D, double *DU, double *B,
		   int *LDB, int *INFO );

#define SINC(A) sin(M_PI * 2.0 * W * (A))/(M_PI * 2.0 * W * (A))
#define NTHREADS 1

mfft::mfft(int _nfft, int _npoints, int _ntapers, const dmatrix &_tapers, const dvector &_lambdas) :
	nfft(_nfft), npoints(_npoints), ntapers(_ntapers)
{}
	
