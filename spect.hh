#ifndef SPECT_H
#define SPECT_H

#include "common.hh"
#include <fftw3.h>
#include "blitz_io.hh"
#include "mtm.hh"

using namespace blitz;

typedef std::complex<double> cmplx;

template <typename T> inline
dmatrix mtmspec(const Array<T,1> &signal, int nfft, double nw, int ntapers,
		int shift, int adapt=1) 
{
	mfft *mtmh = mtm_init_dpss(nfft, nw, ntapers);
	dmatrix out(nfft / 2 + 1, signal.size() / shift, blitz::ColumnMajorArray<2>());
	mtm_spec(mtmh, out.data(), signal.data(), signal.size(), shift, adapt);
	mtm_destroy(mtmh);
	return out;
}

template <typename T> inline
void hamming(Array<T,1> &w)
{
	int n = w.size();
	w = 0.54 - 0.46 * cos (2 * M_PI * tensor::i / (n-1));
}

template <typename T> inline
void hanning(Array<T,1> &w)
{
	int n = w.size();
	w = 0.54 * (1 - cos (2 * M_PI * tensor::i / (n-1)));
}

template <typename T> inline
void gauss(Array<T,1> &w, double sigma)
{
	double n = (double)w.size();
	w = exp(-(sqr((tensor::i - n/2)/sigma))/2.0);
}

template <typename T> inline
void gauss2d(Array<T,2> &w, double sigma1, double sigma2)
{
	double nr = (double)w.rows();
	double nc = (double)w.cols();
	w = exp(-sqr((tensor::i - nr/2)/sigma1)/2.0 - sqr((tensor::j - nc/2)/sigma2)/2.0);
}

template <typename T, int N_rank>
T frac_thresh(const Array<T,N_rank> &data, double frac)
{
	unsigned int i;
	std::vector<T> val(data.size());
	std::copy(data.begin(), data.end(), val.begin());
 	std::sort(val.begin(), val.end());
	double sum = 0, max = (double)blitz::sum(data);
	for (i = 0; (i < val.size()) && (sum < frac); i++)
		sum += (double)val[i] / max;
	return val[i];
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

/*==================================================================
  Complex spectrogram transformation.  Given S(n), n>=0, "smoothing
  window vector" W(0), ... W(N-1), a set of indices [grid], and the
  "within-window location" [r_loc], the spectrogram of S registered on
  [grid] is a matrix, such that its i-th column vector is the Fourier
  coefficients of Z(0),... Z(N-1), where Z(n)=S(x+n-r_loc)*W(n) and
  x=grid(i).  The elements in [grid] can be negative or larger than
  N=total length of S, and for n<0 or n>=N, S(n) is defined to be 0.
  ================================================================*/

template <typename T_in, typename T_win, typename T_out>
void window_signal(const Array<T_in,1> &signal, const Array<T_win,1> &window, 
		   const Array<int,1> &grid, Array<T_out,2> &output, int r_loc=0) {
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


template <typename T1, typename T2>
void mask_spectrogram(const Array<T1,2> &spec, const Array<T2,2> &mask, Array<T1,2> &out)
{
	int n = spec.rows();
	Range all = Range::all();
	out.resize(spec.shape());
	
	out(Range(0,n/2),all) = spec(Range(0,n/2),all) * mask;
	out(Range(n-n/2,n-1),all) = spec(Range(n-n/2,n-1),all) * mask(Range(mask.rows()-1,1,-1),all);
}

class STFT {
private:
	cmatrix buffer;
	fftw_plan fwd,rev;
public:
	STFT(int nfft, int nframes, int nthreads=1);
	~STFT();
	
	template <typename T_in, typename T_win>
	const cmatrix& specgram(const Array<T_in,1> &signal, const Array<T_win,1> &window, 
				const Array<int,1> &grid, int r_loc=0, bool forward=true) 
	{
		window_signal(signal, window, grid, buffer, r_loc);
		execute_transform(forward);
		return buffer;
	}

	template <typename T_in, typename T_win, typename T_out>
	void ispecgram(const Array<T_in,2> &spec, const Array<T_win,1> &window,
		       const Array<int,1> &grid, Array<T_out,1> &output) 
        {
		std::copy(spec.begin(), spec.end(), buffer.begin());
		execute_transform(false);
		overlap_add(window, grid, output);
	}

	template <typename T_win, typename T_out>
	void overlap_add(const Array<T_win,1> &window, const Array<int,1> &grid, Array<T_out,1> &output)
	{
		int N = grid(grid.size()-1) + buffer.rows();
		int Nw = window.size();
		output.resize(N);
		Array<T_win,1> diag(N);

		diag = 0;
		output = 0;
		for (int col = 0; col < buffer.cols(); col++) {
			output(Range(grid(col),grid(col)+Nw)) = window * buffer(Range::all(),col);
			diag(Range(grid(col),grid(col)+Nw)) += blitz::sqr(window);
		}
		output /= diag;
	}

	const cmatrix& get_buffer() const { return buffer; }


protected:
	void execute_transform(bool forward=true);
};

bool is_hermitian(const cvector &vec, double tol=1e-10);

#endif
