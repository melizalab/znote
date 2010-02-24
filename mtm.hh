#ifndef MTM_H
#define MTM_H
/*
 * mtm.hh - C++ library for calculating multitaper spectrograms.
 */
#include <fftw3.h>

/**
 * Multi-taper FFT transformation structure. Contains basic FFT
 * parameters as well as pointers to tapers (i.e. windowing
 * functions), output buffers, and the FFTW plan
 */
typedef struct {
	int nfft;
	int npoints;
	int ntapers;
	double *tapers;
	double *lambdas;
	double *buf;
	fftw_plan plan;
} mfft;

/* initialization and destruction functions */

/**
 * Initialize a multitaper mtm transform using preallocated tapers (i.e. with dpss()))
 *
 * Inputs:
 *   nfft - number of points in the transform
 *   npoints - number of points in the tapers (windows)
 *   ntapers - number of tapers
 *   *tapers - pointer to npoints*ntapers array of windowing functions
 *   *lambdas - eigenvalues for tapers; if NULL, assign weight of 1.0 to each taper
 *
 * Returns:
 *   pointer to mfft_params structure (owned by caller)
 *
 * Note:
 *   pointers to tapers and lambdas are now owned by the return mtfft structure
 */
mfft* mtm_init(int nfft, int npoints, int ntapers, double* tapers, double *lambdas);

/**
 * Initialize a mtfft transform using DPSS tapers (i.e. for a standard
 * multitaper transform)
 *
 * Inputs:
 *   nfft - number of points in the transform/dpss tapers
 *   nw   - time-frequency parameter
 *   ntapers - number of tapers to keep
 *
 * Returns:
 *   pointer to mfft structure (owned by caller)
 */
mfft* mtm_init_dpss(int nfft, double nw, int ntapers);

/**
 * Frees up the mtftt_params structure and dependent data. Note that
 * references to the tapers are considered to be owned by the
 * structure, so if they were calculated elsewhere do not attempt to
 * access them after calling this function.
 */
void mtm_destroy(mfft *mtm);

/* transformation functions */

/**
 * Compute multitaper FFT of a signal. Note that this can be used for
 * single taper FFTs, if the mfft structure has been
 * initialized with a single window
 *
 * Inputs:
 *    mtm - parameters for the transform
 *    data - input data (double-precision floating points)
 *    nbins - the number of time points in the signal
 *
 * Returns:
 *    total power in signal (used in computing adaptive power spectra)
 */
template <typename T>
double mtfft(mfft *mtm, const T *data, int nbins) 
{
	// copy data * tapers to buffer
	int nfft = mtm->nfft;
	int size = mtm->npoints;
	int i,j;
	int nt = (nbins < size) ? nbins : size;
	double pow = 0.0;

	//printf("Windowing data (%d points, %d tapers)\n", nt, mtm->ntapers);
	for (i = 0; i < mtm->ntapers; i++) {
		for (j = 0; j < nt; j++) {
			mtm->buf[j+i*nfft] = mtm->tapers[j+i*size] * data[j];
			pow += data[j] * data[j];
		}
	}

	pow /= mtm->ntapers;
	// zero-pad rest of buffer
	//printf("Zero-pad buffer with %d points\n", mtm->nfft - nt);
	for (i = 0; i < mtm->ntapers; i++) {
		for (j = nt; j < mtm->nfft; j++)
			mtm->buf[j+i*nfft] = 0.0;
	}

	fftw_execute(mtm->plan);

	return pow / nt;
}

/* spectrogram functions */

/**
 * Compute power spectrum from multiple taper spectrogram.  The
 * 'high-res' method is simply the average of the estimates for each
 * taper weighted by the eigenvalue of the taper.  The 'adaptive'
 * method attempts to fit the contribution from each taper to match
 * the total power in the signal.
 *
 * Inputs:
 *   mtm - mfft structure after running mtfft
 *   sigpow - total power in the signal. If zero or less, uses high-res method
 *
 * Outputs:
 *   pow - power spectral density (linear scale) of the signal. Needs to be
 *         preallocated, with dimensions at least nfft/2 + 1;
 */
void mtpower(const mfft *mtm, double *pow, double sigpower);

/**
 *  Compute a multitaper spectrogram by stepping through a signal.
 *  This function 'fills' a spectrogram by calculating the PSD for each
 *  frame in the signal.
 *
 * Inputs:
 *  mtm - mfft structure
 *  samples - input signal
 *  nsamples - number of points in input buffer
 *  shift    - number of samples to shift in each frame
 *  adapt    - if true, use adaptive averaging between tapers (otherwise 'high-res')
 *
 * Outputs:
 *  spec     - reassigned spectrogram. needs to be allocated and zero-filled before calling
 *
 */
template <typename T>
void mtm_spec(mfft *mtm, double *spec, const T *samples, int nsamples, int shift, int adapt) 
{
	int t;
	int nbins = nsamples / shift;
	int real_count = mtm->nfft / 2 + 1;
	double sigpow;

	for (t = 0; t < nbins; t++) {
		sigpow = mtfft(mtm, samples+(t*shift), nsamples-(t*shift));
		mtpower(mtm, spec+(t*real_count), (adapt) ? sigpow : 0.0);
	}
}


#endif
