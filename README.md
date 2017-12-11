# znote

**znote** consists of a suite of libraries and programs for
spectrotemporal decomposition of acoustic signals. Many bioacoustic
signals, such as birdsong, consist of a set of spectrotemporally
disjoint "objects". These objects may overlap in time with each other,
as when a bird is vocalizing with both sides of its syrinx
simultaneously; or overlap with sounds from other animals or the
environment, as in field recordings. It is often possible to separate
these objects from each other in the spectrotemporal domain, but it then
becomes necessary to reconstruct the sound pressure waveforms giving
rise to the isolated spectrotemporal components.

The first stage is to identify the components of the vocalization.
**znote~label~** provides one method of doing this, by finding all the
connected components in a signal. Briefly, the spectrogram of the signal
is computed, and all the points which are above this spectrogram are
grouped into contiguous features. The output is an array with the same
dimensions of the spectrogram, in which each time-frequency point is
given an integer code indicating which component it belongs to.

The second stage is to invert the spectrographic transform for each of
the identified components. This is done by **znote~extract~**, which
takes as input the original signal and the array indicating which points
in the spectrogram belong to which features. It is not necessary to use
**znote~label~** to generate this array, and the feature-file can be
manipulated after it is generated to group or split components.
**zedit**, a MATLAB GUI for simple feature manipulation is provided.

For more information on the algorithm, see:

Meliza CD, Chi Z, Margoliash D, "Representations of Conspecific Song by
Starling Secondary Forebrain Auditory Neurons: Towards a Hierarchical
Framework". J Neurophysiology, [doi:10.1152/jn.00464.2009](http://doi.org/10.1152/jn.00464.2009)

License and warranty
--------------------

Use of the code is free for non-commercial purposes under the Creative
Commons Attribution-Noncommercial-Share Alike 3.0 United States License
(<http://creativecommons.org/licenses/by-nc-sa/3.0/us/>).

C Daniel Meliza, Zhiyi Chi, and Daniel Margoliash (or the above paper)
will be acknowledged as the source of the algorithms in any publications
reporting its use or the use of any modified version of the program.

THE PROGRAMS ARE PROVIDED "AS IS" WITHOUT WARRANTY OF MERCANTABILITY OR
FITNESS FOR A PARTICULAR PURPOSE OR ANY OTHER WARRANTY, EXPRESS OR
IMPLIED. IN NO EVENT SHALL THE UNIVERSITY OF CHICAGO, THE UNIVERSITY OF
CONNECTICUT, OR DRS. MELIZA, CHI OR MARGOLIASH BE LIABLE FOR ANY DIRECT
OR CONSEQUENTIAL DAMAGES RESULTING FROM USE OF THE PROGRAMS. THE USER
BEARS THE ENTIRE RISK FOR USE OF THE PROGRAMS.

To contact the authors, please consult the official repository for znote
at <http://www.github.org/melizalab/znote>

Compilation and Installation
----------------------------

### Dependencies

**znote~label~** and **znote~extract~** require

-   a modern C++ compiler
-   scons &gt;= 1.3
-   libsndfile &gt;= 1.0.17 (<http://www.mega-nerd.com/libsndfile>)
-   blitz &gt;= 0.9 (<http://blitz.sourceforge.net/>)
-   fftw &gt;= 3.3 (<http://www.fftw.org>)
-   LAPACK (specifically functions dsterf and dgtsv)

\*zedit\* requires a reasonably recent version of MATLAB

### Mac OS X

Znote was developed on OS X 10.6 (Snow Leopard). OS X comes with LAPACK
as part of the vecLib framework. FFTW, blitz, and libsndfile can be
installed through the macports system, as follows:

``` bash
sudo port install blitz fftw-3 libsndfile
```

Scons is also available through macports, although this will cause a
separate copy of python to be installed under the macports tree, if it
is not already.

``` bash
sudo port install scons
```

The SConstruct file, as supplied, should compile the programs on any
system with XCode installed:

``` bash
scons -Q -j2
```

### Linux

Znote has been compiled and tested on an x86~64~ Opteron system running
Redhat EL5, with the ATLAS BLAS libraries and NETLIB LAPACK libraries
installed under `/usr/local`. If necessary, edit the lib~path~ and
include~path~ in SConstruct to point to directories where the
appropriate headers and libraries are located.

### Windows

Znote has been successfully compiled on windows XP using the MinGW
compiler (i.e. gcc 3.4). Getting all the dependencies compiled is
tricky. Blitz (0.9) does not compile, or work with gcc 4.4. I recommend
installing the enthought Python distribution, which ships with MinGW,
and installing Msys or Cygwin in order to run the configure scripts for
the dependencies. FFTW3, blitz and libsndfile should compile out of the
box. Use static libraries if possible to avoid having dependencies on
DLLs. LAPACK needs to be compiled with g77, not gfortran. The SConstruct
file assumes that the dependencies have been installed under the znote
source tree.

Hopefully, you will not have to do any of this - znote~label~.exe and
znote~extract~.exe were compiled on Windows XP. Like the rest of znote,
no support is offered.

### Performance notes

FFTW can be configured to use multiple threads, which can offer some
speed improvement for large sound files. For small sound files, the
overhead of setting up multiple threads is rarely worth it. To enable
multiple threads in znote, edit SConstruct and set threads to some
number less than or equal to the number of cores in your system.

Usage
-----

An example is provided under `test/`, consisting of a recorded starling
vocalization (`A8.wav`) and labeled components of the vocalizations
(`A8_feats.bin`). To test your install, run

``` {.example}
znote_label --pad --recon test/A8.wav test/A8_feats.bin
```

This will generate a separate file for each component of the signal, and
a reconstruction formed by adding all the components together. Compare
`A8_recon.wav` with `A8_recon_reference.wav`.

### znote~label~

``` {.example}
znote_label [--nfft <i>] [--fftshift <i>]  [--ntapers <i>] [--nw <f>]
                [--thresh <f>] [--df <f>] [--dt <f>]
                [--min-size <f>] <input>
```

`nfft`: controls the size of the FFT analysis window. Default 512, which
is appropriate most signals sampled at around 44 kHz. Larger values give
higher frequency resolution at the expense of lower temporal resolution.
The value of nfft is most important at this stage, because it determines
the time-frequency resolution of algorithm that detects connected
components.

`fftshift`: controls the spacing between FFT analysis windows. Default
is 10, which gives a substantial amount of overlap between frames.
Increasing the value can increase the speed of the algorithm, at some
cost to the temporal resolution during labelling.

`nw`: this program uses a multitaper algorithm to estimate spectral
density. Increasing the time-bandwidth product increases thes stability
of these estimates, but at the expense of lower spectral resolution. The
default value of 3.5 gives a decent amount of smoothing. Larger values
give more smoothing, but neighboring components may get smeared
together. Smaller values can improve resolution between neighboring
components, but tend to underestimate the ST extent of the components
and increase the number of points where the power goes above threshold
spuriously. Needs to be a half-integer (i.e. 3,3.5,4,...)

`ntapers`: provides further control over spectrogram estimation.
Defaults to nw/2-1, which is generally considered to be the optimal
value.

`thresh`: set the minimum power for a component. This can be specified
in absolute terms, in dB, or relative to the total amount of power in
the signal. If the value is greater than 1.0, the threshold is
calculated as an absolute value, and only the points in the spectrogram
where the power is greater than this value are considered to be "above
water" for the detection of components. If less than 1.0, the absolute
threshold is calculated as the power corresponding to the quantile
&lt;thresh&gt;. Default is 0.5 (or 50%). Note that the relative
threshold is calculated on a linear scale, so 50% of the power is often
confined to a fairly small portion of the signal.

`df`: control frequency resolution of component search algorithm.
Components are considered to be connected if they are less than df Hz
apart. Defaults to 200 Hz. Along with dt, increasing values lead to
fewer, larger components.

`dt`: control temporal resolution of component search algorithm.
Defaults to 2 ms.

`min-size`: Components with less than &lt;min-size&gt; kHz-ms area are
dropped.

The input file to **znote~label~** can be a sound file (in any format
libsndfile understands), or a .bin file containing the spectrogram of
the signal. Consult blitz~io~.hh for documentation on the .bin format.
The behaviors of many of the flags change when using a pre-calculated
spectrogram, so this is not recommended for novice users.

The program outputs a .bin file indicating which points in the
spectrogram belong to which features.

### znote~extract~

``` {.example}
znote_extract [--fbdw <f>] [--tbdw <f>]
               [--feat <i>] [--pad] [--del] [--recon]
               SIGNAL LABELFILE
```

**znote~extract~** uses the labels defined in `LABELFILE` to generate
masks, which it uses to extract the associated time series in `SIGNAL`.
The masks are generated with a Gaussian roll-off filter, the parameters
of which are controlled on the command line:

`fbdw`: Set frequency bandwidth for Gaussian roll-off mask. Defaults to
200 Hz. Larger values reduce edge effects, but at the cost of
potentially interfering with neighboring components, or including more
noise.

`tbdw`: Set time bandwidth for smoothing kernel. Defaults to 2 ms.

`feat`: By default, the program extracts all the component defined in
&lt;labels&gt;; set this value to a nonnegative integer to restrict to a
single component.

`pad`: By default, the program generates unpadded output files; if this
flag is set, then the output signals are the same length as the input
signal, with all points where the component was not present set to 0.

`del`: If set, the program will also generate deletions, which are
calculated by substracting (at the appropriate temporal offset) the
extracted components from the original signal.

`recon`: If set, the program will sum all the extracted components at
their original offset and output the resulting sum.

`SIGNAL` must be a sound file, because the program needs the original
phase information to reconstruct the signals.

`LABELFILE` can be any integer bin file, including the file output by
**znote~label~**. The dimensions of the file will be used to control the
FFT parameters of the extraction algorithm.

Output:

**znote~extract~** writes one wave file for each extracted component. If
the input file is named signal.wav, the output files will be named
signal~feature000~.wav, signal~feature001~.wav, etc.

For component deletions, the output files are named as
signal~fdel000~.wav, etc

The reconstruction has the name signal~recon~.wav

### zedit

**zedit** is a simple MATLAB interface for editing .bin files. It allows
merging and splitting of components while visualizing the spectrogram of
the corresponding signal. To edit components for a signal, run zedit in
MATLAB as follows:

``` {.example}
>> zedit <wavefile>
```

zedit runs **znote~label~** to generate spectrograms and calculate
connected components. If the executable is not in your path, you may
need to edit zedit~params~.m When the program first runs, it will
calculate the spectrogram of &lt;wavefile&gt; and display it with a
single contour indicating where the threshold lies.

The parameters of the spectrographic transform can be changed in the
FFT/MTM panel. The threshold value can be edited manually or by clicking
on the colorbar to the right of the spectrogram.

In the LabelSet panel, to calculate components, click the Label button.
Note: this will overwrite the file &lt;wavefile&gt;~labels~.bin. To load
a previously generated label file, click Load.

When a labelset is selected, a list of features is displayed in the
Features panel. Selecting one or more features causes them to be
displayed in the spectrogram. Features can be merged with the Merge
button, or split by clicking the Lasso button. After clicking Lasso,
click points on the spectrogram to define a polygon around the feature
of interest. Click the middle mouse button to close the polygon and
split the feature. Only currently selected features are affected.

Save the edited labelset by clicking Save in the LabelSet panel. Choose
a name for the output file; this can be used with **znote~extract~** to
generate the signals associated with the components.

### Best practices

The algorithm in **znote~extract~** can generate artifacts near edges of
signals. If this is a problem, be sure to pad your signals with silence
on either end.

Another source of error is when components overlap. Overlap is caused by
the rolloff filter, so if your components are too close together, try
decreasing these values. You can tell if overlap is occurring when
**znote~extract~** outputs a line like "Max feature overlap: 1.90596".
If this value is 1.0 there is no overlap.

Version History
---------------

### 1.1.0

First public release.

### 1.2.0

Updated to compile with blitz 0.10. Should still compile with blitz 0.9.
Fixed a bug where the mask was not being correctly applied to the row
corresponding to half Nyquist.
