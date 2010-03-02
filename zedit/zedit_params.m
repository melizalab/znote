function S = zedit_params()
%
% ZEDIT_PARAMS returns a structure with default zedit parameters
%
% Copyright C Daniel Meliza, Z Chi 2010.  Licensed for use under Creative
% Commons Attribution-Noncommercial-Share Alike 3.0 United States
% License (http://creativecommons.org/licenses/by-nc-sa/3.0/us/).


S = struct('nw', 3.5, 'nfft', 512, 'fftshift', 10, ...
           'df', 200, 'dt', 2, 'thresh', 60, 'min_size', 4000, 'fbdw', ...
           200, 'tbdw', 2, 'labeller','znote_label');

% note: edit last entry in S to set path to znote_label program
