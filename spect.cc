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

void
STFT::overlap_add(const dvector &window, const ivector &grid, dvector &output,
                  int first_col, int last_col)
{
        if (last_col < first_col) last_col = grid.size();
        int Nw = window.size();
        int N = grid(last_col-1) - grid(first_col) + Nw;
        output.resize(N);
        dvector diag(N);
        // std::cout << window << std::endl;

        diag = 0;
        output = 0;
        for (int col = first_col; col < last_col; col++) {
                int t = grid(col) - grid(first_col);
                output(Range(t,t+Nw-1)) += window * blitz::real(buffer(Range::all(),col));
                // std::cout << output(Range(t,t+Nw-1)) << std::endl;
                diag(Range(t,t+Nw-1)) += blitz::sqr(window);
        }
        output /= blitz::where(diag>0,diag,1.0);
}

void
STFT::ispecgram()
{
        fftw_execute(rev);
        buffer /= buffer.rows();
}

void
STFT::ispecgram(const cmatrix &spec)
{
        assert(spec.rows()==buffer.rows());
        assert(spec.cols()==buffer.cols());
        std::copy(spec.begin(), spec.end(), buffer.begin());
        ispecgram();
}


double
STFT::ispecgram(const cmatrix &spec, const dmatrix &mask)
{
        assert(spec.rows()==buffer.rows());
        assert(spec.cols()==buffer.cols());
        assert(mask.rows()==spec.rows()/2+1);
        int n = spec.rows();
        double maxpow = 0;
        for (int c = 0; c < buffer.cols(); c++) {
                buffer(0,c) = spec(0,c) * mask(0,c);
                for (int r = 1; r < mask.rows()-1; r++) {
                        buffer(r,c) = spec(r,c) * mask(r,c);
                        buffer(n-r,c) = spec(n-r,c) * mask(r,c);
                        maxpow = fmax(maxpow, abs(buffer(r,c)));
                }
                if (n % 2 == 0) {
                        buffer(n/2,c) = spec(n/2,c) * mask(mask.rows()-1,c);
                }
        }
        ispecgram();
        return maxpow;
}

