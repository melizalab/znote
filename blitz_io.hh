#ifndef BLITZIO_H
#define BLITZIO_H
/**
 * @file   blitz_io.hh
 * @author Daniel Meliza <dmeliza@uchicago.edu>
 * @date   Mon Mar  1 13:35:27 2010
 * 
 * @brief  C++ functions for reading and writing data in Blitz++ arrays.
 * 
 * Copyright C Daniel Meliza, Z Chi 2010.  Licensed for use under Creative
 * Commons Attribution-Noncommercial-Share Alike 3.0 United States
 * License (http://creativecommons.org/licenses/by-nc-sa/3.0/us/).
 */
#include <stdexcept>
#include <cstdio>
#include <blitz/array.h>
#include <sndfile.hh>

/** 
 * Read .bin file, which is a binary storage format for 1D or 2D array
 * data.
 * 
 * @param filename - input file
 * @param T - type of data in the file
 * @param out - output array (contents are overwritten)
 * 
 * @return - number of samples read
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

/** 
 * Write array to a .bin file. The first two integers in the file
 * indicate the data dimensions, and the remainder of the file
 * consists of the contents of the array, written in whatever storage
 * order the array is using.
 * 
 * @param filename - output file name
 * @param data - output data array; can be of any type; behavior undefined for rank>2
 * 
 * @return - number of samples written
 */
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
	

/**
 * A 1D array with an associated sample rate. Defines some functions
 * for reading and writing to sound files via libsndfile.
 */
template<typename T>
struct timeseries {

	blitz::Array<T,1> samples; /**< Array holding the samples */
	unsigned int samplerate; /**< The sampling rate of the time series */

/** 
 * Construct a new timeseries from an existing 1D array
 * 
 * @param T Type of the data
 * @param _samples The sample values. This array is *referenced* by the timeseries object
 * @param _samplerate The sampling rate of the time series
 */	
	timeseries(blitz::Array<T,1> &_samples, unsigned int _samplerate) :
		samples(_samples), samplerate(_samplerate) {}

/** 
 * Construct a new, empty timeseries
 * 
 * @param _nsamples Number of samples
 * @param _samplerate Sampling rate
 */
	timeseries(unsigned int _nsamples, unsigned int _samplerate) :
		samples(blitz::Array<T,1>(_nsamples)), samplerate(_samplerate) {}

/** 
 * Construct a timeseries by reading in a sound file
 * 
 * @param filename The name of the file to read
 */
	timeseries(const char *filename) { read_pcm(filename);}
	timeseries(const std::string &filename) { read_pcm(filename.c_str());}

/** 
 * Read data from a sound file into the object.  The reference to the
 * current data is lost and replaced with a new array, which the
 * object owns.
 * 
 * @param filename The name of the sound file
 * 
 * @return the number of samples read
 */
	int read_pcm(const char *filename) {
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
		return count;
	}

/** 
 * Write the array to a sound file.
 * 
 * @param filename Name of the output file
 * @param format Format of the output file - default is 16-bit little-endian WAV
 * @param norm If true (default), normalize floating point data [-1.0,1.0] becomes the range of the output format.
 * 
 * @return the number of sample written
 */
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
		
