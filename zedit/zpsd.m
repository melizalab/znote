function [PSD,T,F] = zpsd(signalfile, varargin)
%
% ZPSD Computes the spectrotemporal PSD of an input waveform
% using multitaper methods.
%
% [PSD,T,F] = ZPSD(signalfile, [parameters])
%
% Returns the ST power spectrum in a matrix.
%
% Copyright C Daniel Meliza, Z Chi 2010.  Licensed for use under Creative
% Commons Attribution-Noncommercial-Share Alike 3.0 United States
% License (http://creativecommons.org/licenses/by-nc-sa/3.0/us/).

error(nargchk(1,2,nargin));
params = mergestruct(zedit_params, varargin{:});
PSD = [];
[signal,Fs,bd] = wavread(signalfile);

cstr = sprintf('LD_LIBRARY_PATH="" %s --nfft %d --fftshift %d --nw %1.1f --df %3.2f --dt %3.2f --thresh %3.2f --min-size %3.2f --spec %s', ...
               params.labeller, params.nfft, params.fftshift, ...
               params.nw, params.df, params.dt, ...
               params.thresh, params.min_size, signalfile);

fprintf('Computing labels: %s\n', cstr);
[status, result] = system(cstr);
if status==0
  fprintf('Done!\n');
  % now load the binfiles
  [p f e] = fileparts(signalfile);
  specfile = sprintf('%s_spec.bin', f);
  PSD   = bimatrix(specfile, 'double','F');
  PSD   = log10(max(PSD,1.0)) * 10;
  T = linspace(0,Fs/2,size(PSD,1));
  F = linspace(0,size(PSD,2)*params.fftshift*1000/Fs, ...
             size(PSD,2));
else
  fprintf('Error running label command: %s\n', result);
end
