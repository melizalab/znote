#include <complex>
#include "spect.hh"

STFT::STFT(int nfft, int nframes, int nthreads)
{
	TinyVector<int,2> buf_shape(nfft, nframes);
	cmplx *out_data = (cmplx*)fftw_malloc(nfft * nframes * sizeof(cmplx));
	buffer.reference(cmatrix(out_data, buf_shape, neverDeleteData));
	fftw_complex *buf = reinterpret_cast<fftw_complex*>(buffer.data());
	ivector n(nframes);
	n = nfft;

	//fftw_plan_with_nthreads(nthreads);
	fwd  = fftw_plan_many_dft(1, n.data(), nframes,
				  buf, NULL, nframes, 1,
				  buf, NULL, nframes, 1,
				  FFTW_FORWARD, FFTW_MEASURE);
	rev  = fftw_plan_many_dft(1, n.data(), nframes,
				  buf, NULL, nframes, 1,
				  buf, NULL, nframes, 1,
				  FFTW_BACKWARD, FFTW_MEASURE);
}

STFT::~STFT()
{
	fftw_destroy_plan(fwd);
	fftw_destroy_plan(rev);
	fftw_free(buffer.data());
}

void
STFT::execute_transform(bool forward)
{
	fftw_execute(forward ? fwd : rev);
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
