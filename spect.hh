#ifndef SPECT_H
#define SPECT_H

#include "math.hh"
#include <fftw3.h>
#include "mtm.hh"

using namespace blitz;

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
Array<T,1> hamming(unsigned int n)
{
	Array<T,1> out(n);
	out = (0.54 - 0.46 * cos (2 * M_PI * tensor::i / (n-1)));
	return out;
}

template <typename T> inline
Array<T,1> hanning(unsigned int n)
{
	Array<T,1> out(n);
	out = 0.54 * (1 - cos (2 * M_PI * tensor::i / (n-1)));
	return out;
}

template <typename T, int N_rank>
T frac_thresh(const Array<T,N_rank> &data, double frac)
{
	int i;
	std::vector<T> val(data.size());
	std::copy(data.begin(), data.end(), val.begin());
 	std::sort(val.begin(), val.end());
	double sum = 0, max = (double)blitz::sum(data);
	for (i = 0; (i < val.size()) && (sum < frac); i++)
		sum += (double)val[i] / max;
	return val[i];
}
	

// template <typename T>
// Array<std::complex<double>,2> specgram(const Array<T,1> &signal, 

#endif
