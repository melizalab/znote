/**
 * @file   spect.cc
 * @author Daniel Meliza <dmeliza@uchicago.edu>
 * @date   Mon Mar  1 13:38:31 2010
 * 
 * Copyright C Daniel Meliza, Z Chi 2010.  Licensed for use under Creative
 * Commons Attribution-Noncommercial-Share Alike 3.0 United States
 * License (http://creativecommons.org/licenses/by-nc-sa/3.0/us/).
 * 
 */
#include "spect.hh"
#include <complex>

#if THREADS>1
static int t = fftw_init_threads();
#endif


STFT::STFT(int nfft, int nframes, int nthreads)
{
	TinyVector<int,2> buf_shape(nfft, nframes);
	cmplx *out_data = (cmplx*)fftw_malloc(nfft * nframes * sizeof(cmplx));
	buffer.reference(cmatrix(out_data, buf_shape, neverDeleteData));
	fftw_complex *buf = reinterpret_cast<fftw_complex*>(buffer.data());
	ivector n(nframes);
	n = nfft;

	int planner_flag = (nframes > 5000) ? FFTW_MEASURE : FFTW_ESTIMATE;

#if THREADS>1
	fftw_plan_with_nthreads(nthreads);
#endif
	fwd  = fftw_plan_many_dft(1, n.data(), nframes,
				  buf, NULL, nframes, 1,
				  buf, NULL, nframes, 1,
				  FFTW_FORWARD, planner_flag);  // change to FFTW_MEASURE for huge files
	rev  = fftw_plan_many_dft(1, n.data(), nframes,
				  buf, NULL, nframes, 1,
				  buf, NULL, nframes, 1,
				  FFTW_BACKWARD, planner_flag);
}

STFT::~STFT()
{
	fftw_destroy_plan(fwd);
	fftw_destroy_plan(rev);
	fftw_free(buffer.data());
}

bool
is_hermitian(const cvector &vec, double tol)
{
	int n = vec.size();
	
	if (fabs(imag(vec(0))) > tol) return false;
	if ((n % 2)==0 && fabs(imag(vec(n/2))) > tol) return false;
	cvector pos(vec(Range(1,n/2)));
	cvector neg(vec(Range(n-1,n-n/2,-1)));
	if (blitz::any(blitz::abs(real(pos)-real(neg))>tol)) return false;
	if (blitz::any(blitz::abs(imag(pos)+imag(neg))>tol)) return false;
	return true;
}
