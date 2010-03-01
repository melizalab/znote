/**
 * @file   mtm.cc
 * @author Daniel Meliza <dmeliza@uchicago.edu>
 * @date   Mon Mar  1 13:38:31 2010
 * 
 * Copyright C Daniel Meliza, Z Chi 2010.  Licensed for use under Creative
 * Commons Attribution-Noncommercial-Share Alike 3.0 United States
 * License (http://creativecommons.org/licenses/by-nc-sa/3.0/us/).
 * 
 */
#include <cstdio>
#include <cmath>
#include <float.h>
#include <complex>
extern "C" {
#include "dpss.h"
}
#include "mtm.hh"

mfft*
mtm_init(int nfft, int npoints, int ntapers, double* tapers, double *lambdas)
{
	mfft *mtm;
	int *n_array, i;
	fftw_r2r_kind *kind;
	mtm = (mfft*)malloc(sizeof(mfft));

	mtm->nfft = nfft;
	mtm->npoints = npoints;
	mtm->ntapers = ntapers;
	mtm->tapers = tapers;
	if (lambdas)
		mtm->lambdas = lambdas;
	else {
		mtm->lambdas = (double*)malloc(ntapers*sizeof(double));
		for (i = 0; i < ntapers; i++) mtm->lambdas[i] = 1.0;
	}

	mtm->buf = (double*)fftw_malloc(nfft*ntapers*sizeof(double));

	// set up fftw plan
	n_array = new int[ntapers];
	kind = new fftw_r2r_kind[ntapers];
	for (i = 0; i < ntapers; i++) {
		n_array[i] = nfft;
		kind[i] = FFTW_R2HC;
	}

	mtm->plan = fftw_plan_many_r2r(1, n_array, ntapers,
				       mtm->buf, NULL, 1, nfft,
				       mtm->buf, NULL, 1, nfft,
				       kind, FFTW_MEASURE);

	delete n_array;
	delete kind;
	return mtm;
}


void
mtm_destroy(mfft *mtm)
{
	if (mtm->plan) fftw_destroy_plan(mtm->plan);
	if (mtm->tapers) free(mtm->tapers);
	if (mtm->lambdas) free(mtm->lambdas);
	if (mtm->buf) fftw_free(mtm->buf);
	free(mtm);
}


void
mtpower(const mfft *mtm, double *pow, double sigpow)
{
	int nfft = mtm->nfft;
	int ntapers = mtm->ntapers;
	int real_count = nfft / 2 + 1;
	int imag_count = (nfft+1) / 2;  // not actually the count but the last index
	int t,n;

	if (sigpow<=0.0 || ntapers==1) {
		memset(pow, 0, real_count*sizeof(double));
		for (t = 0; t < ntapers; t++) {
			for (n = 0; n < real_count; n++)
				pow[n] += mtm->buf[t*nfft+n]*mtm->buf[t*nfft+n]*mtm->lambdas[t]/ntapers;
			for (n = 1; n < imag_count; n++) {
				pow[n] += mtm->buf[t*nfft+(nfft-n)]*mtm->buf[t*nfft+(nfft-n)]*mtm->lambdas[t]/ntapers;
			}
		}
	}
	else {
		double est, num, den, w;
		double tol, err;
		double *Sk;
		Sk = (double*)calloc(ntapers*real_count, sizeof(double));
		for (t = 0; t < ntapers; t++) {
			for (n = 0; n < real_count; n++)
				Sk[t*real_count+n] += mtm->buf[t*nfft+n]*mtm->buf[t*nfft+n]*mtm->lambdas[t];
			for (n = 1; n < imag_count; n++)
				Sk[t*real_count+n] += mtm->buf[t*nfft+(nfft-n)]*mtm->buf[t*nfft+(nfft-n)]*mtm->lambdas[t];
			//Sk[t*nfft+n] *= 2;
		}
		// initial guess is average of first two tapers
		err = 0;
		for (n = 0; n < real_count; n++) {
			pow[n] = (Sk[n] + Sk[real_count+n])/2;
			err += abs(pow[n]);
		}

		tol = 0.0005 * sigpow / nfft;
		err /= nfft;
		//printf("err: %3.4g; tol: %3.4g\n", err, tol);
		//for(t = 0; t < ntapers; t++)
		//      printf("%3.4g ", sigpow * (1 - mtm->lambdas[t]));
		//printf("\n");
		while (err > tol) {
			err = 0;
			for (n = 0; n < real_count; n++) {
				est = pow[n];
				num = den = 0;
				//printf("%d: est=%3.4g; ", n, est);
				for (t=0; t < ntapers; t++) {
					w = est / (est * mtm->lambdas[t] + sigpow * (1 - mtm->lambdas[t]));
					w = w * w * mtm->lambdas[t];
					//printf("%3.4g ",Sk[t*real_count+n]);
					num += w * Sk[t*real_count+n];
					den += w;
				}
				pow[n] = num/den;
				err += fabs(num/den-est);
			}
			//printf("err: %3.4g\n", err);
		}
		free(Sk);
	}
	// adjust power for one-sided spectrum
	for (n = 1; n < imag_count; n++)
		pow[n] *= 2;
}



mfft*
mtm_init_dpss(int nfft, double nw, int ntapers)
{
	double *tapers, *lambdas;
	tapers = (double*)malloc(nfft*ntapers*sizeof(double));
	lambdas = (double*)malloc(nfft*sizeof(double));
	dpss(tapers, lambdas, nfft, nw, ntapers);
	return mtm_init(nfft, nfft, ntapers, tapers, lambdas);
}
