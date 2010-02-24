
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <complex.h>
#include <fftw3.h>
#include "dpss.h"


/* some LAPACK prototypes */
extern void dsterf_(int *N, double *D, double *E, int *INFO);

extern void dgtsv_(int *N, int *NRHS,
		   double *DL, double *D, double *DU, double *B,
		   int *LDB, int *INFO );

#define SINC(A) sin(M_PI * 2.0 * W * (A))/(M_PI * 2.0 * W * (A))


/**
 * Scale a vector by its L2 norm
 *
 * Inputs:
 *   N - number of points
 *   x - vector; altered in place
 */
static void
renormalize(int N, double *x)
{
	int i;
	double norm = 0.0;
	for (i = 0; i < N; i++)
		norm += x[i]*x[i];

	norm = sqrt(norm);
	for (i = 0; i < N; i++)
		x[i] /= norm;
}


/**
 * Compute the self-convolution of a vector using FFT.
 *
 * Inputs:
 *   N - number of points
 *   x - input vector
 *
 * Outputs:
 *   y - output vector (not allocated; needs to have N points)
 */
static void
fftconv(int N, const double *x, double *y)
{
	int i;
	double *X;
	fftw_complex *X1, *X2;
	fftw_plan plan;

	//fftw_init_threads();

	X = (double*)calloc(N * 2, sizeof(double));
	X1 = (fftw_complex*)fftw_malloc((N+1) * sizeof(fftw_complex));
	X2 = (fftw_complex*)fftw_malloc((N+1) * sizeof(fftw_complex));

	memcpy(X, x, sizeof(double)*N);
	plan = fftw_plan_dft_r2c_1d(N*2, X, X1, FFTW_ESTIMATE);
	fftw_execute(plan);

	// flip x and compute again
	for (i = 0; i < N; i++)
		X[i] = x[N-1-i];
	plan = fftw_plan_dft_r2c_1d(N*2, X, X2, FFTW_ESTIMATE);
	fftw_execute(plan);

	for (i = 0; i < N; i++)
		X1[i] *= X2[i];

	// inverse fft
	plan = fftw_plan_dft_c2r_1d(N*2, X1, X, FFTW_ESTIMATE);
	fftw_execute(plan);

	for (i = 0; i < N; i++)
		y[i] = X[i] / N / 2;
	fftw_free(X1);
	fftw_free(X2);
	free(X);
	fftw_destroy_plan(plan);
}

int
dpss(double *tapers, double *lambda, int npoints, double NW, int k)
{
	int i, j, m, rv;
	double *d, *sd, *dd1, *dd2, *ee1, *ee2;
	double *taper;

	double W, ff;

	if ((NW < 0) || (k < 1) || (k >= npoints) || (npoints < 0) || (NW >= npoints/2))
		return -1;

	W = NW/npoints;

	d = (double*)malloc(npoints*sizeof(double));
	sd = (double*)malloc(npoints*sizeof(double));
	dd1 = (double*)malloc(npoints*sizeof(double));
	dd2 = (double*)malloc(npoints*sizeof(double));
	ee1 = (double*)malloc((npoints)*sizeof(double));
	ee2 = (double*)malloc((npoints)*sizeof(double));

	for (i = 0; i < npoints; i++) {
		ff = (npoints - 1 - 2*i);
		d[i] = dd1[i] = 0.25 * cos(2*M_PI*W) * ff * ff;
		sd[i] = ee1[i] = (i+1) * (npoints-(i+1))/2.0;
	}

	// lapack eigenvalue solver; values stored in d in increasing order
	dsterf_(&npoints,d,ee1,&rv);
	if (rv != 0) return -2;

	// set up tridiagonal equations:
	for (j = 0; j < k; j++) {
		taper = tapers + j * npoints;  // point into tapers array
		lambda[j] = d[npoints-(j+1)];
		// initialize taper
		for (i = 0; i < npoints; i++)
			taper[i] = sin((j+1) * M_PI * i / (npoints-1));

		for (m = 0; m < 3; m++) {
			// all inputs destroyed by dgtsv
			for (i = 0; i < npoints; i++) {
				dd2[i] = dd1[i] - lambda[j];
				ee1[i] = ee2[i] = sd[i];
			}
			i = 1;
			dgtsv_(&npoints, &i, ee1, dd2, ee2, taper, &npoints, &rv);
			if (rv != 0) return -2;
			renormalize(npoints, taper);
		}

		// fix sign of taper
		if ((j+1) % 2==1) {
			// calculate sum
			ff = 0.0;
			for (i = 0; i < npoints; i++)
				ff += taper[i];
			if (ff < 0.0) {
				for (i = 0; i < npoints; i++)
					taper[i] *= -1;
			}
		}
		else if (taper[2] < 0.0) {
			for (i = 0; i < npoints; i++)
				taper[i] *= -1;
		}

		// calculate lambdas
		fftconv(npoints, taper, dd2);

		ff = 2.0 * W * dd2[npoints-1];  // last point
		for (i = 0; i < npoints-1; i++)
			ff += dd2[i] * 4.0 * W * SINC(npoints-1-i);

		lambda[j] = ff;

	}
	free(d);
	free(sd);
	free(dd1);
	free(dd2);
	free(ee1);
	free(ee2);
	return 0;
}

