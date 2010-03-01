/**
 * @file   znote_label.cc
 * @author Daniel Meliza <dmeliza@uchicago.edu>
 * @date   Mon Mar  1 13:38:31 2010
 * 
 * Copyright C Daniel Meliza, Z Chi 2010.  Licensed for use under Creative
 * Commons Attribution-Noncommercial-Share Alike 3.0 United States
 * License (http://creativecommons.org/licenses/by-nc-sa/3.0/us/).
 * 
 */
#include "common.hh"
#include "blitz_io.hh"
#include "mtm.hh"
#include "components.hh"
#include <getopt.h>
#include <iostream>

using namespace std;

static const string program_name = "znote_label";
static const string program_version = "1.1.0";

double spec_threshold = 0.5;
int nfft = 512;
int fft_shift = 10;
int ntapers = 5;
double nw = 3.5;

double f_nbhd_r = 200.0;
double t_nbhd_r = 2.0;
double min_area = 4000.0;

string file_name;

void
usage()
{
	cout << program_name << " (version " << program_version << ")\n\n ";
	cout << program_name << " [--nfft <i>] [--fftshift <i>]  [--ntapers <i>] [--nw <f>]\n"
	     << "                 [--thresh <f>] [--df <f>] [--dt <f>]\n"
	     << "                 [--min-size <f>] <input>\n" << endl;
	     
  
	cout << "<i> indicates an integer argument; <f> a float. See documentation for details\n"
	     << "<input> can be either a sound file or a spectrogram.\n"
	     << "Output is to <pcmfile>_labels.bin" << endl;

	exit(-1);
}


void
parse_args(int argc, char **argv)
{
	for (int i = 1; i < argc; i++) {
		if (strncmp(argv[i], "--nfft",6)==0)
			nfft = atoi(argv[++i]);
		else if (strncmp(argv[i], "--fftshift",10)==0)
			fft_shift = atoi(argv[++i]);
		else if (strncmp(argv[i], "--ntapers",9)==0)
			ntapers = atoi(argv[++i]);
		else if (strncmp(argv[i], "--nw",4)==0)
			nw = atof(argv[++i]);
		else if (strncmp(argv[i], "--df", 4) == 0)
			f_nbhd_r = atof(argv[++i]);
		else if (strncmp(argv[i], "--dt", 4) == 0)
			t_nbhd_r = atof(argv[++i]);
		else if (strncmp(argv[i], "--min-size", 10) == 0)
			min_area = atof(argv[++i]);
		else if (strncmp(argv[i], "--thresh", 8) == 0)
			spec_threshold = atof(argv[++i]);
		else
			file_name = string(argv[i]);
	}
	       
	if (file_name.size() < 1) {
		cerr << "Error: must supply a file for input..." << endl;
		exit(-1);
	}
}

int
main(int argc, char **argv) {

	if (argc == 1)
		usage();
	parse_args(argc, argv);

	dmatrix spec;
 	coord_list nbhd;
	imatrix labels;
	int freq_r, time_r, min_size;
	string froot, fext;

	cout << "* Program: " << program_name << endl
	     << "* Version: " << program_version << endl
	     << "* Input: " << file_name << endl;
	
	if (splitext(file_name, froot, fext) < 0) {
		cout << "* ERROR: Unable to determine input file type" << endl;
		exit(-1);
	}
	if (fext!=".bin") {
		timeseries<short> pcm(file_name.c_str());
		cout <<  "* Samples: " << pcm.samples.size() << endl
		     <<  "* Samplerate: " << pcm.samplerate << endl
		     <<  "* Nfft: " << nfft << " samples" << endl
		     <<  "* Shift: " << fft_shift << " samples" << endl
		     <<  "* Tapers: " << ntapers << endl
		     <<  "* Time-freq product: " << nw << endl
		     <<  "* Minimum feature area: " << min_area * 0.001 << " Hz-ms " << endl;
	
		spec.reference(mtmspec(pcm.samples, nfft, nw, ntapers, fft_shift));
		freq_r = (int)(f_nbhd_r * nfft / pcm.samplerate);
		time_r = (int)(t_nbhd_r * pcm.samplerate / fft_shift * .001);
		min_size = (int)(min_area * nfft / fft_shift * .001);
		
	}
	else {
		read_bin(file_name.c_str(), spec);
		freq_r = (int)(f_nbhd_r);
		time_r = (int)(t_nbhd_r);
		min_size = (int)min_area;
	}
	cout << "* Spectrogram dimensions: [" << spec.rows() << " " << spec.cols() << "]" << endl;

	imatrix tspec(spec.shape());
	if (spec_threshold >= 1.0) {
		cout << "* Feature threshold: " << spec_threshold << endl;
		tspec = (blitz::log10(spec) > (spec_threshold / 10));
	}
	else {
		cout << "* Relative threshold: " << spec_threshold << endl;
		double abs_thresh = frac_thresh(spec, spec_threshold);
		cout << "* Absolute threshold: " << log10(abs_thresh) * 10 << endl;
		tspec = (spec > abs_thresh);
	}

	cout << "* Frequency search radius: " << freq_r << endl
	     << "* Time search radius: " << time_r << endl;

	ellipse_neighborhood(nbhd, freq_r, time_r);

 	clist_vector clist(matrix_cbranches(tspec, nbhd, labels));
	cout << "* Extracted features: " << clist.size() << endl
	     << "* Minimum size: " << min_size << endl 
	     << "------------------------------------" << endl
	     << "feat\tpixels\trow\tcol" << endl;
	clist_vector::const_iterator it = clist.begin();
	for (int i = 0; it != clist.end(); it++, i++) {
		blitz::RectDomain<2> fbounds = component_bounds(*it);
		cout << i << "\t" << it->size() << "\t";
		cout << fbounds.lbound(0) << "\t";
		cout << fbounds.lbound(1) << "\t";
		if (it->size() < min_size) {
			delete_component(labels, *it);
			cout << "dropped";
		}
		cout << endl;
	}
	renumber_components(labels);
	cout << "------------------------------------" << endl
	     << "* Total features: " << blitz::max(labels) << endl;

 	write_bin((froot + "_labels.bin").c_str(),labels);
	
}

