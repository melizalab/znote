function [PSD,T,F] = fedit_psd(pcmfile, varargin)
%
% FEDIT_PSD Computes the spectrotemporal PSD of an input waveform
% using multitaper methods.
%
% [PSD,T,F] = FEDIT_PSD(pcmfile, [parameters])
%
% Returns the ST power spectrum in a matrix, with columns
% for each time point. If no output argument is supplied,
% plots the PSD in the current figure.
%
% CDM, 8/2006

error(nargchk(1,2,nargin));
params = mergestruct(fog_params, varargin{:});

cstr = sprintf('%s -t 100 --nfft %d --fftshift %d --mtm %f %s',...
               params.fog_psd, params.nfft, params.fftshift, params.mtm, ...
               pcmfile);
fprintf('Generating DB spectro: %s\n', cstr);
[status, result] = system(cstr);
fprintf('Done!\n');

% Load the spectrogram
[p f e] = fileparts(pcmfile);
dbfile = sprintf('%s_mtm_%d_%d_%.1f.bin', f, params.nfft, params.fftshift, ...
                 params.mtm);
PSD = bimatrix(dbfile, 'double');
% the values < 0 really aren't worth much and they eff up the colormap
PSD(PSD < 0) = 0;
T = linspace(0,params.Fs/2,size(PSD,1));
F = linspace(0,size(PSD,2)*params.fftshift*1000/params.Fs, ...
             size(PSD,2));

% if the function is being called without return value, plot
% the data
if nargout == 0
     imagesc(F,T,PSD*10);
     xlabel('Time (ms)');ylabel('F (Hz)');
     axis xy, axis tight, box on;
     %colorbar
end
