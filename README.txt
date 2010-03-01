
znote - identify and extract connected spectrotemporal components

Znote consists of a suite of libraries and programs for
spectrotemporal decomposition of acoustic signals.  Many bioacoustic
signals, such as birdsong, consist of a set of spectrotemporally
disjoint "objects".  These objects may overlap in time with each
other, as when a bird is vocalizing with both sides of its syrinx
simultaneously; or overlap with sounds from other animals or the
environment, as in field recordings.  It is often possible to separate
these objects from each other in the spectrotemporal domain, but it
then becomes necessary to reconstruct the sound pressure waveforms
giving rise to the isolated spectrotemporal components.

The first stage is to identify the components of the vocalization.
znote_label provides one method of doing this, by finding all the
connected components in a signal.  Briefly, the spectrogram of the
signal is computed, and all the points which are above this
spectrogram are grouped into contiguous features.  The output is an
array with the same dimensions of the spectrogram, in which each
time-frequency point is given an integer code indicating which
component it belongs to.

The second stage is to invert the spectrographic transform for each of
the identified components.  This is done by znote_extract, which takes
as input the original signal and the array indicating which points in
the spectrogram belong to which features.  It is not necessary to use
znote_label to generate this array, and the feature-file can be
manipulated after it is generated to group or split components.  A
MATLAB GUI for simple feature manipulation is provided.

For more information on the algorithm, see:

[1] Meliza CD, Chi Z, Margoliash D, "Representations of Conspecific
Song by Starling Secondary Forebrain Auditory Neurons: Towards a
Hierarchical Framework". J Neurophysiology, doi:10.1152/jn.00464.2009

Use of the code is licensed for non-commercial purposes under the
Creative Commons Attribution-Noncommercial-Share Alike 3.0 United
States License (http://creativecommons.org/licenses/by-nc-sa/3.0/us/). 
For other purposes, contact the authors at:

* Dan Meliza (dan@meliza.org)
* Zhiyi Chi (zchi@merlot.stat.uconn.edu)
* Dan Margoliash (d-margoliash@uchicago.edu)

******** Installation

znote depends on the Blitz C++ array library for data manipulation,
the FFTW library for computing FFT transforms, and the libsndfile
library for reading acoustic signal data.  It also requires the LAPACK
fucntions dsterf and dgtsv to generate the discrete prolate spherical
sequences used by znote_label
