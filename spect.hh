#ifndef SPECT_H
#define SPECT_H
/**
 * @file   spect.hh
 * @author Daniel Meliza <dmeliza@uchicago.edu>
 * @date   Mon Mar  1 13:33:47 2010
 *
 * @brief  Functions to compute forward and inverse spectrograms.
 *
 * Copyright C Daniel Meliza, Z Chi 2010.  Licensed for use under Creative
 * Commons Attribution-Noncommercial-Share Alike 3.0 United States
 * License (http://creativecommons.org/licenses/by-nc-sa/3.0/us/).
 */
#include "common.hh"
#include <fftw3.h>
#include <algorithm>
#include "blitz_io.hh"

using namespace blitz;

typedef std::complex<double> cmplx;

template <typename T> inline
void hamming(Array<T,1> &w)
{
	firstIndex i;
	int n = w.size();
	w = 0.54 - 0.46 * cos (2 * M_PI * i / (n-1));
}

template <typename T> inline
void hanning(Array<T,1> &w)
{
	firstIndex i;
	int n = w.size();
	w = 0.54 * (1 - cos (2 * M_PI * i / (n-1)));
}

template <typename T> inline
void gauss(Array<T,1> &w, double sigma)
{
	firstIndex i;
	w = exp(-(sqr((i - w.size()/2)/sigma))/2.0);
}

template <typename T> inline
void gauss2d(Array<T,2> &w, double sigma1, double sigma2)
{
	firstIndex i;
	secondIndex j;
	w = exp(-sqr((i - w.rows()/2)/sigma1)/2.0 - sqr((j - w.cols()/2)/sigma2)/2.0);
}

inline
int freq2row(double f_hbdw, int nfft, int sample_rate)
{
	return (int)(f_hbdw * nfft / sample_rate);
}

inline
double row2freq(int rows, int nfft, int sample_rate)
{
	return ((double)rows * sample_rate / nfft);
}

inline
int time2col(double t_hbdw, int fft_shift, int sample_rate)
{
	return (int)(t_hbdw * sample_rate / fft_shift * .001);
}

inline
double col2time(int cols, int fft_shift, int sample_rate)
{
	return (double)cols * fft_shift / sample_rate * 1000;
}

/**
 * Calculate windowed signal in multiple frames.  Prepares a signal
 * for transformation into a spectrogram by sliding an analysis window
 * along the signal and multiplying it by a windowing function. Given
 * S(n), n>=0, "smoothing window vector" W(0), ... W(N-1), a set of
 * indices [grid], and the "within-window location" [r_loc], the
 * spectrogram of S registered on [grid] is a matrix, such that its
 * i-th column vector is the Fourier coefficients of Z(0),... Z(N-1),
 * where Z(n)=S(x+n-r_loc)*W(n) and x=grid(i).  The elements in [grid]
 * can be negative or larger than N=total length of S, and for n<0 or
 * n>=N, S(n) is defined to be 0.
 *
 * @param signal - input signal (S)
 * @param window - smoothing window (W)
 * @param grid - grid of analysis time points
 * @param output - output array (Z)
 * @param r_loc - relative location in the window
 */
template <typename T_in, typename T_win, typename T_out>
void window_signal(const Array<T_in,1> &signal, const Array<T_win,1> &window,
		   const Array<int,1> &grid, Array<T_out,2> &output, int r_loc=0) {
#if BZ_MAJOR_VERSION > 0 || BZ_MINOR_VERSION > 9
        using blitz::extrema::min;
#endif
	int nframes = min(grid.size(), output.cols());
	int nfreq = min(window.size(), output.rows());
	int i,c,lp,rp;
	for (c = 0; c < nframes; c++) {
		lp = grid(c) - r_loc;
		rp = lp + nfreq;
		for (i = lp; i < rp; i++) {
			if (i < 0 || i >= signal.size())
				output(i-lp,c) = 0;
			else
				output(i-lp,c) = signal(i) * window(i-lp);
		}
	}
}


/**
 * Wrapper around FFTW for computing spectrograms.
 */
class STFT {
private:
	cmatrix buffer;		/**< internal buffer for forward and backward transforms */
	fftw_plan fwd,rev;
public:

/**
 * Initialize STFT object for transforms of a given size and number.
 *
 * @param nfft Number of points in the transform
 * @param nframes Number of frames to transform
 * @param nthreads Number of threads to use in transform
 */	STFT(int nfft, int nframes, int nthreads=1);
	~STFT();

/**
 * Compute spectrogram of a signal.
 *
 * @param signal Input signal
 * @param window Smoothing window
 * @param grid Grid points for each frame
 * @param r_loc Relative location of window in each frame
 *
 * @return reference to internal buffer after transformation
 */	template <typename T_in, typename T_win>
	const cmatrix& specgram(const Array<T_in,1> &signal, const Array<T_win,1> &window,
				const Array<int,1> &grid, int r_loc=0)
	{
		window_signal(signal, window, grid, buffer, r_loc);
		fftw_execute(fwd);
		return buffer;
	}

/**
 * Compute inverse transform of data in the buffer.
 */
	void ispecgram();

/**
 * Compute inverse spectrogram.
 *
 * @param spec 2D complex array - copied into internal buffer
 */
	void ispecgram(const cmatrix &spec);

/**
 * Compute inverse spectrogram of masked data.
 *
 * @param spec 2D complex array, dimensions (N x L)
 * @param mask mask to apply to array as it's being copied in, dimensions (N/2+1 x L)
 *
 * @return maximum power in the masked array
 */
	double ispecgram(const cmatrix &spec, const dmatrix &mask);

/**
 * Calculate real-valued signal using overlap and add method. Assumes
 * ispecgram() has been called first.  For normal behavior, the
 * smoothing window and time grid should be the same as used to
 * calculate the original spectrogram, but some interesting effects
 * can be achieved, by warping the grid for example.
 *
 * @param window Smoothing window (used to calculate original spectrogram)
 * @param grid Frame time grid (used to calculate original spectrogram)
 * @param output 1D double-valued output. Resized to appropriate size (losing reference to original contents)
 * @param first_col first analysis column (default 0)
 * @param last_col last analysis column (default last)
 */
	void overlap_add(const dvector &window, const ivector &grid, dvector &output,
			 int first_col=0, int last_col=-1);

	const cmatrix& get_buffer() const { return buffer; }

};

/**
 * Calculate whether a complex vector is hermitian.  Given a vector C
 * with N samples, C(0) and C(N/2) [if N is even] will be real-valued,
 * and C(i) = C*(N-i), where * indicates the complex conjugate.
 *
 * @param vec Input vector
 * @param tol Tolerance of computation
 *
 * @return true if hermitian, false otherwise
 */
bool is_hermitian(const cvector &vec, double tol=1e-10);

#endif
