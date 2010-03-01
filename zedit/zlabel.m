function [L,PSD] = fedit_label(pcmfile, varargin)
%
% FEDIT_LABEL  Helper to FEDIT GUI for labelling components in a
% spectrogram.
%
% [L,PSD] = FEDIT_LABEL(pcmfile, params)   
%
% The pcmfile is passed to fog_label with the following parameters
% from the params struct:
%
% .nfft:      The number of samples in the FFT/MTM window
% .fftshift:  The amount to shift each window (controls time
%             resolution of spectrogram)
%
% .mtm:       The time-bandwidth resolution of the MTM process.
% .thresh:    The threshhold, in dB, used to convert the
%             spectrogram into a binary matrix.
%
% .df:        The radius, in the frequency dimension, of the
%             neighborhood function used for component labelling.
%             Units of Hz.
%
% .dt:        Like dt, but in the time dimension (units of ms)
%
% TODO:  allow dynamic specification of the neighborhood function
% (needs an update to fog_label)

error(nargchk(1,2,nargin));
params = mergestruct(fog_params, varargin{:});

% by setting the threshhold unreasonably high we get no features
cstr = sprintf('%s -v -t %3.2f -s %d --nfft %d --fftshift %d --mtm %1.1f --df %d --dt %d %s', ...
               params.fog_label, params.thresh, params.min_size, params.nfft, params.fftshift, ...
               params.mtm, params.df, params.dt, ...
               pcmfile);

fprintf('Computing labels: %s\n', cstr);
[status, result] = system(cstr);
fprintf('Done!\n');

% now load the binfiles
[p f e] = fileparts(pcmfile);
dbfile = sprintf('%s_mtm_%d_%d_%.1f.bin', f, params.nfft, params.fftshift, ...
                 params.mtm);
lblfile = sprintf('%s_labels.bin', f);

PSD = bimatrix(dbfile, 'double');
PSD(PSD < 0) = 0;
L   = bimatrix(lblfile, 'int') + 1;
