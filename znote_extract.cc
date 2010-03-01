
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

	char buf[128];
	imatrix labels;
	string sfroot, sfext;
	int nfft, fft_shift, pad_col;

	cout << "* Program: " << program_name << endl
	     << "* Version: " << program_version << endl
	     << "* Signal: " << signal_file << endl
	     << "* Labels: " << label_file << endl;

	if (splitext(signal_file, sfroot, sfext) < 0) {
		cout << "* ERROR: Unable to determine signal file type" << endl;
		exit(-1);
	}

	timeseries<double> pcm(signal_file);
	cout << "* Samples: " << pcm.samples.size() << endl
	     << "* Samplerate: " << pcm.samplerate << endl;
	
	read_bin(label_file, labels);
	nfft = (labels.rows() - labels.rows() % 2) * 2;
	fft_shift = int(pcm.samples.size() / labels.cols());
	clist_vector features(get_components(labels));
	cout << "* Label array dimensions: [" << labels.rows() << " " << labels.cols() << "]" << endl
	     << "* Features: " << features.size() << endl
	     << "* Nfft: " << nfft << endl
	     << "* Shift: " << fft_shift << endl;
	pad_col = (int)(120.0 / fft_shift);

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
	cmatrix spec(stft.get_buffer().copy());

	// Gaussian roll-off filter
	int fn = 1 + 2 * freq2row(f_hbdw, nfft, pcm.samplerate);
	int tn = 1 + 2 * time2col(t_hbdw, fft_shift, pcm.samplerate);
	cout << "* Frequency rolloff: " << f_hbdw << " Hz (" << fn << " bins)" << endl
	     << "* Time rolloff: " << t_hbdw << " ms (" << tn << " bins)" << endl;
	dmatrix gfilt(fn, tn);
	gauss2d(gfilt, double(fn)/4, double(tn)/4);
	cout << gfilt << endl;
	
	dmatrix masked_tot(labels.shape());
	mask_sum(features,gfilt,coord(fn/2,tn/2),masked_tot);
	masked_tot = blitz::max(masked_tot,1.0);
	
	ivector feat_nums(1);
	if (feat_num > -1)
		feat_nums(0) = feat_num;
	else
		arange(feat_nums,0,features.size());

	cout << "-----------------------------------------" << endl
	     << "feat" << '\t' << "t.onset" << '\t' << "f.onset" << '\t' 
	     << "t.size" << '\t' << "f.size" << '\t'
	     << "maxDB" << '\t' << "area" << '\t' << "samples" << endl;

	dmatrix mask(labels.shape());
	dvector recon(pcm.samples.size());
	for (int i = 0; i < feat_nums.size(); i++) {

		cout << feat_nums(i);
		coord_list feature = features[feat_nums(i)];

		RectDomain<2> fbounds = component_bounds(feature);
		cout << '\t' << fbounds.lbound(1) //col2time(tmin,fft_shift,pcm.samplerate)
		     << '\t' << fbounds.lbound(0) //row2freq(fmin,nfft,pcm.samplerate)
		     << '\t' << fbounds.ubound(1) - fbounds.lbound(1) 
		     << '\t' << fbounds.ubound(0) - fbounds.lbound(0);
		
 		mask = 0;
 		make_mask(feature,gfilt,coord(fn/2,tn/2),mask);
 		mask /= masked_tot;

 		cmatrix masked;
 		mask_spectrogram(spec, mask, masked);
		cout << '\t' << log10(blitz::max(abs(masked))) * 10
		     << '\t' << feature.size();

 		dvector output;
		int start_col, stop_col;
		if (pad_features) {
			start_col = 0;
			stop_col = grid.size()-1;
		}
		else {
			start_col = max(0,fbounds.lbound(1)-pad_col);
			stop_col = min(grid.size()-1,fbounds.ubound(1)+pad_col);
		}
		stft.ispecgram(masked, window, grid, output, start_col, stop_col);

		sprintf(buf, "%s_feature_%03d.wav", sfroot.c_str(), feat_nums(i));
 		
		cout << '\t' << timeseries<double>(output, pcm.samplerate).write_pcm(buf);
 		if (output_deletions) {
 			dvector deletion = pcm.samples.copy();
			deletion(Range(grid(start_col),grid(stop_col))) -= output;
 			sprintf(buf, "%s_fdel_%03d.wav", sfroot.c_str(), feat_nums(i));
			timeseries<double>(deletion, pcm.samplerate).write_pcm(buf);
 		}
		if (output_recon)
			recon(Range(grid(start_col),grid(stop_col))) += output;
			
 		cout << endl;
	}
	if (output_recon) {
		sprintf(buf, "%s_recon.wav", sfroot.c_str());
		cout << "* Wrote reconstruction to: " << buf << endl;
		timeseries<double>(recon, pcm.samplerate).write_pcm(buf);
	}
}

