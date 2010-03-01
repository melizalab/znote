#ifndef BLITZIO_H
#define BLITZIO_H

#include <stdexcept>
#include <cstdio>
#include <blitz/array.h>
#include <sndfile.hh>

/* 
 * Define some functions for reading and writing data in blitz arrays
 */

template<typename T>
int read_bin(const char *filename, blitz::Array<T,2> &out) {
	FILE *fp;
	int count, expected;
	int shape[2];

	if ((fp = fopen(filename, "rb")) == NULL)
		throw std::runtime_error("File does not exist");
	fread(shape, sizeof(int), 2, fp);
	blitz::TinyVector<int,2> tshape(shape[0], shape[1]);
	expected = shape[0]*shape[1];

	out.resize(tshape);
	count = fread(out.data(), sizeof(T), expected, fp);
	fclose(fp);
	return count;
}

template<typename T> inline
int read_bin(const std::string &filename, blitz::Array<T,2> &out) {
	return read_bin(filename.c_str(), out);
}


template<typename T, int N>
int write_bin(const char *filename, const blitz::Array<T,N> &data) {
	FILE *fp;
	int nr, nc, ct;
	const T *pData = data.data();

	if ((fp = fopen(filename, "wb")) == NULL)
		return 0;
	nr = data.rows();
	nc = data.cols();
	fwrite(&nr, sizeof(int), 1, fp);
	fwrite(&nc, sizeof(int), 1, fp);
	ct = fwrite(pData, sizeof(T), nr*nc, fp);
	fclose(fp);
	return ct;
}


template<typename T, int N> inline
int write_bin(const std::string &filename, const blitz::Array<T,N> &data) {
	return write_bin(filename.c_str(), data);
}
	

template<typename T>
struct timeseries {

	blitz::Array<T,1> samples;
	unsigned int samplerate;

	timeseries(blitz::Array<T,1> &_samples, unsigned int _samplerate) :
		samples(_samples), samplerate(_samplerate) {}

	timeseries(unsigned int _nsamples, unsigned int _samplerate) :
		samples(blitz::Array<T,1>(_nsamples)), samplerate(_samplerate) {}

	timeseries(const char *filename) {
		sf_count_t count, nframes;
		SndfileHandle fp = SndfileHandle(filename);
		if (fp.error())
			throw std::runtime_error("Sound file does not exist");
		if (fp.channels() != 1)
			throw std::runtime_error("Sound file has more than 1 channel");
		samplerate = fp.samplerate();
		nframes = fp.frames();
		T *data = new T[nframes];
		count = fp.read(data, nframes);
		samples.reference(blitz::Array<T,1>(data, nframes, blitz::deleteDataWhenDone));
		samples.resize(count);
	}

	timeseries(const std::string &filename) {
		timeseries(filename.c_str());
	}

	int write_pcm(const char *filename, int format = SF_FORMAT_WAV | SF_FORMAT_PCM_16, 
		      bool norm=true) {
		SndfileHandle fp = SndfileHandle(filename, SFM_WRITE, format, 1, samplerate);
		if (fp.error())
			throw std::runtime_error("Error opening file: " + std::string(fp.strError()));
		fp.command(SFC_SET_NORM_DOUBLE,NULL,norm ? SF_TRUE : SF_FALSE);
		fp.command(SFC_SET_NORM_FLOAT,NULL,norm ? SF_TRUE : SF_FALSE);
		const T *pData = samples.data();
		return fp.writef(pData, samples.size());
	}
};

#endif
		
