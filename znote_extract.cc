
#include "common.hh"
#include "blitz_io.hh"
#include "spect.hh"
#include "components.hh"
#include <getopt.h>
#include <iostream>

using namespace std;

static const string program_name = "znote_extract";
static const string program_version = "1.1.0";

// The half bandwidth of the smoothing kernel in the frequency domain, in HZ
double f_hbdw = 200.0;

// The half bandwidth of the smoothing kernel in the temporal domain; in msec
double t_hbdw = 2.0;

// the user can choose to only extract one feature. If negative, extract all
int feat_num = -1;

bool output_deletions = false;
bool output_recon = false;
bool pad_features = false;

string signal_file;
string label_file;

void
usage()
{
	cout << program_name << " (version " << program_version << ")\n\n ";
	cout << program_name << " [--fbdw <f>] [--tbdw <f>]\n"
	     << "                 [--feat <i>] [--pad] [--del] [--recon]\n"
	     << "                 <signal> <labels>\n" << endl;
	     
  
	cout << "<i> indicates an integer argument; <f> a float. See documentation for details\n"
	     << "<signal> must be either a sound file.\n"
	     << "The dimensions of <labels> control spectrographic calculations" << endl;


	exit(-1);
}


void
parse_args(int argc, char **argv)
{
	for (int i = 1; i < argc; i++) {
		if (strncmp(argv[i], "--fbdw",6)==0)
			f_hbdw = atof(argv[++i]);
		else if (strncmp(argv[i], "--tbdw",6)==0)
			t_hbdw = atof(argv[++i]);
		else if (strncmp(argv[i], "--feat",6)==0)
			feat_num = atoi(argv[++i]);
		else if (strncmp(argv[i], "--del",5)==0)
			output_deletions = true;
		else if (strncmp(argv[i], "--pad",5)==0)
			pad_features = true;
		else if (strncmp(argv[i], "--recon", 7) == 0)
			output_recon = true;
		else if (signal_file.size()==0)
			signal_file.append(argv[i]);
		else if (label_file.size()==0)
			label_file.append(argv[i]);
	}
	       
	if (signal_file.size() < 1) {
		cerr << "Error: must supply a signal file for input..." << endl;
		exit(-1);
	}
	if (label_file.size() < 1) {
		cerr << "Error: must supply a label file for input..." << endl;
		exit(-1);
	}
}


int
main(int argc, char **argv) {

	if (argc == 1)
		usage();
	parse_args(argc, argv);

	dmatrix spec;
	imatrix labels;
	string sfroot, sfext;
	int nfft, fft_shift;

	cout << "* Program: " << program_name << endl
	     << "* Version: " << program_version << endl
	     << "* Signal: " << signal_file << endl
	     << "* Labels: " << label_file << endl;

	if (splitext(signal_file, sfroot, sfext) < 0) {
		cout << "* ERROR: Unable to determine signal file type" << endl;
		exit(-1);
	}

	timeseries<short> pcm(signal_file.c_str());
	cout << "* Samples: " << pcm.samples.size() << endl
	     << "* Samplerate: " << pcm.samplerate << endl;
	
	read_bin(label_file.c_str(), labels);
	nfft = (labels.rows() - labels.rows() % 2) * 2;
	fft_shift = int(pcm.samples.size() / labels.cols());
	clist_vector features(get_components(labels));
	cout << "* Label array dimensions: [" << labels.rows() << " " << labels.cols() << "]" << endl
	     << "* Features: " << features.size() << endl
	     << "* Nfft: " << nfft << endl
	     << "* Shift: " << fft_shift << endl;

	if (features.size()==0) {
		cout << "* Error: No valid features defined in " << label_file << endl;
		exit(-1);
	}
	if (feat_num > (int)features.size()) {
		cout << "* Error: " << label_file << " does not define feature " << feat_num << endl;
		exit(-1);
	}

	dvector window(nfft);
	hanning(window);

	ivector grid;
	arange(grid,0,pcm.samples.size(),fft_shift);

 	STFT stft(nfft, grid.size());
	stft.specgram(pcm.samples, window, grid);
	write_bin((sfroot + "_spec.bin").c_str(),stft.get_buffer());

	// Gaussian roll-off filter
	int fn = 1 + 2 * freq2row(f_hbdw, nfft, pcm.samplerate);
	int tn = 1 + 2 * time2col(t_hbdw, fft_shift, pcm.samplerate);
	cout << "* Frequency rolloff: " << f_hbdw << " Hz (" << fn << " bins)" << endl
	     << "* Time rolloff: " << t_hbdw << " ms (" << tn << " bins)" << endl;
	dmatrix gfilt(fn, tn);
	gauss2d(gfilt, double(fn)/4, double(tn)/4);
	
	dmatrix masked_tot(labels.shape());
	mask_sum(features,gfilt,coord(fn/2,tn/2),masked_tot);
	masked_tot = blitz::max(masked_tot,1.0);
	
	ivector feat_nums(1);
	if (feat_num > 0)
		feat_nums(0) = feat_num;
	else
		arange(feat_nums,0,features.size());

	cout << "-----------------------------------------" << endl
	     << "feat" << '\t' << "t.onset" << '\t' << "f.onset" << '\t' 
	     << "t.size" << '\t' << "f.size" << '\t'
	     << "maxDB" << '\t' << "area" << endl;

	dmatrix mask(labels.shape());
	for (int i = 0; i < feat_nums.size(); i++) {

		char buf[128];
		cout << feat_nums(i);
		coord_list feature = features[feat_nums(i)];

		RectDomain<2> fbounds = component_bounds(feature);
		cout << '\t' << fbounds.lbound(0) //col2time(tmin,fft_shift,pcm.samplerate)
		     << '\t' << fbounds.lbound(1) //row2freq(fmin,nfft,pcm.samplerate)
		     << '\t' << fbounds.ubound(0) - fbounds.lbound(0) 
		     << '\t' << fbounds.ubound(1) - fbounds.lbound(1);
		
 		mask = 0;
 		make_mask(feature,gfilt,coord(fn/2,tn/2),mask);
 		mask /= masked_tot;
		

 		cmatrix masked;
 		mask_spectrogram(stft.get_buffer(), mask, masked);
 		sprintf(buf, "%s_feature_%03d.bin", sfroot.c_str(), feat_nums(i));
 		write_bin(buf, masked);
		cout << '\t' << log10(blitz::max(abs(masked))) * 10
		     << '\t' << feature.size() << endl;

// 		dvector output;
// 		//stft.ispecgram(masked, window, grid, output);
 		cout << endl;
	}
}

