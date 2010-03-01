#ifndef MTM_H
#define MTM_H
/**
 * @file   mtm.hh
 * @author Daniel Meliza <dmeliza@uchicago.edu>
 * @date   Mon Mar  1 13:33:47 2010
 * 
 * @brief  C++ library for calculating multitaper spectrograms.
 * 
 * Copyright C Daniel Meliza, Z Chi 2010.  Licensed for use under Creative
 * Commons Attribution-Noncommercial-Share Alike 3.0 United States
 * License (http://creativecommons.org/licenses/by-nc-sa/3.0/us/).
 */
#include <blitz/array.h>
#include <fftw3.h>

/**
 * Multi-taper FFT transformation structure. Contains basic FFT
 * parameters as well as pointers to tapers (i.e. windowing
 * functions), output buffers, and the FFTW plan
 */
typedef struct {
	int nfft;		/**< number of points in the transform */
	int npoints;		/**< number of points in the tapers */
	int ntapers;		/**< number of tapers */
	double *tapers;		/**< pointer to tapers (npoints x ntapers) */
	double *lambdas;	/**< pointer to taper weights (ntapers) */
	double *buf;		/**< transform buffer (npoints x ntapers) */
	fftw_plan plan;		/**< fftw plan */
} mfft;

/* initialization and destruction functions */

/**
 * Initialize a multitaper mtm transform using preallocated tapers
 * (i.e. with dpss())). NB: pointers to tapers and lambdas are now
 * owned by the return mtfft structure
 *
 * @param nfft - number of points in the transform
 * @param npoints - number of points in the tapers (windows)
 * @param ntapers - number of tapers
 * @param *tapers - pointer to npoints*ntapers array of windowing functions
 * @param *lambdas - eigenvalues for tapers; if NULL, assign weight of 1.0 to each taper
 *
 * @return pointer to mfft_params structure (owned by caller)
 *   
 */
mfft* mtm_init(int nfft, int npoints, int ntapers, double* tapers, double *lambdas);

/**
 * Initialize a mtfft transform using DPSS tapers (i.e. for a standard
 * multitaper transform)
 *
 * @param nfft number of points in the transform/dpss tapers
 * @param nw time-frequency parameter (should be a half-integer, e.g. 3.5)
 * @param ntapers number of tapers to keep (should be no more than nw/2-1)
 *
 * @return pointer to mfft structure (owned by caller)
 */
mfft* mtm_init_dpss(int nfft, double nw, int ntapers);

/**
 * Frees up the mtftt_params structure and dependent data. Note that
 * references to the tapers are considered to be owned by the
 * structure, so if they were calculated elsewhere do not attempt to
 * access them after calling this function.
 *
 * @param mtm  structure to free up
 */
void mtm_destroy(mfft *mtm);

/* transformation functions */

/**
 * Compute multitaper FFT of a signal. Note that this can be used for
 * single taper FFTs, if the mfft structure has been
 * initialized with a single window
 *
 * @param mtm    parameters for the transform
 * @param data   input data (double-precision floating points)
 * @param nbins  the number of time points in the signal
 *
 * @return total power in signal (used in computing adaptive power spectra)
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

	for (i = 0; i < mtm->ntapers; i++) {
		for (j = 0; j < nt; j++) {
			mtm->buf[j+i*nfft] = mtm->tapers[j+i*size] * data[j];
			pow += data[j] * data[j];
		}
	}

	pow /= mtm->ntapers;
	// zero-pad rest of buffer
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
 * @param mtm    mfft structure after running mtfft
 * @param sigpow total power in the signal. If zero or less, uses high-res method
 * @param pow    output, power spectral density (linear scale) of the signal. Needs to be
 *               preallocated, with dimensions at least nfft/2 + 1;
 */
void mtpower(const mfft *mtm, double *pow, double sigpower);

/**
 *  Compute a multitaper spectrogram by stepping through a signal.
 *  This function 'fills' a spectrogram by calculating the PSD for each
 *  frame in the signal.
 *
 *  @param mtm      mfft structure
 *  @param samples  input signal
 *  @param nsamples number of points in input buffer
 *  @param shift    number of samples to shift in each frame
 *  @param adapt    if true, use adaptive averaging between tapers (otherwise 'high-res')
 *
 *  @param spec     output, reassigned spectrogram. needs to be allocated and zero-filled before calling
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


/** 
 * Compute multitaper spectrogram of a signal
 * 
 * @param signal 1D array of real-valued samples
 * @param nfft number of points in the analysis window
 * @param nw time-frequency product (should be half-integer)
 * @param ntapers number of tapers (should be no more than nw*2-1)
 * @param shift number of points between analysis frames
 * @param adapt compute adaptive PSDs (1) or high-res (0)
 * 
 * @return array of doubles with dimensions (nfft/2+1, npoints/shift)
 */
template <typename T>
blitz::Array<double,2> mtmspec(const blitz::Array<T,1> &signal, int nfft, double nw, int ntapers,
			       int shift, int adapt=1) 
{
	mfft *mtmh = mtm_init_dpss(nfft, nw, ntapers);
	blitz::Array<double,2> out(nfft / 2 + 1, signal.size() / shift, blitz::ColumnMajorArray<2>());
	mtm_spec(mtmh, out.data(), signal.data(), signal.size(), shift, adapt);
	mtm_destroy(mtmh);
	return out;
}


#endif
