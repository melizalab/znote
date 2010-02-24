
#include "math.hh"
#include "blitz_io.hh"
#include "spect.hh"
#include "components.hh"
#include <iostream>

using namespace std;
using namespace blitz::tensor;

int
main(int argc, char **argv) {

	string infile = "test.wav";

	timeseries<short> pcm(infile.c_str());
	cout << "Read " << pcm.samples.size() << " samples from " << infile << endl;
	
	dmatrix spec = mtmspec(pcm.samples, 512, 3.5, 6, 10);

 	coord_list nbhd;
 	ellipse_neighborhood(nbhd, 5, 5);

	imatrix tspec(blitz::log10(spec) > 6.0);
	imatrix components;
	write_bin("ttest.bin",tspec);
	matrix_cbranches(tspec, nbhd, components);
	write_bin("mtest.bin",components);
	
}
