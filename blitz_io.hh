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
class bin_reader {
private:
	FILE *fp;
public:
	bin_reader(const char *filename) {
		fp = fopen(filename, "rb");
	}
	~bin_reader() {
		fclose(fp);
	}
	blitz::Array<T,2> read() {
		int count, expected;
		int shape[2];
		fread(shape, sizeof(int), 2, fp);
		blitz::TinyVector<int,2> tshape(shape[0], shape[1]);
		expected = shape[0]*shape[1];

		T *data = new T[expected];
		count = fread(data, sizeof(T), expected, fp);

		if (count < expected)
			fprintf(stderr, "Only read %d values, expected %d; wrong data type?",
				count, expected);
		blitz::Array<T,2> out(data, tshape, blitz::deleteDataWhenDone);
		rewind(fp);
		return out;
	}
};

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
	

template<typename T>
struct timeseries {

	blitz::Array<T,1> samples;
	unsigned int samplerate;

	timeseries(blitz::Array<T,1> &_samples, unsigned int _samplerate) :
		samples(_samples), samplerate(_samplerate) {}

	timeseries(const char *filename) {
		sf_count_t count, nframes;
		SndfileHandle fp = SndfileHandle(filename);
		if (fp.channels() != 1)
			throw std::runtime_error("Sound file has more than 1 channel");
		samplerate = fp.samplerate();
		nframes = fp.frames();
		T *data = new T[nframes];
		count = fp.read(data, nframes);
		samples.reference(blitz::Array<T,1>(data, nframes, blitz::deleteDataWhenDone));
		samples.resize(count);
	}

	int write_pcm(const char *filename, int format) {
		SndfileHandle fp = SndfileHandle(filename, SFM_WRITE, format, 1, samplerate);
		const T *pData = samples.data();
		return fp.write(pData, samples.size());
	}
};
	

#endif
		
