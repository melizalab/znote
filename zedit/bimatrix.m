function M=bimatrix(fname, precision,order);
% BIMATRIX Loads binary data from a file
%
% M = BIMATRIX(fname, precision,order)
%
% Read a matrix from a binary file.  The first two numbers in the
% file are of integer*4 type, specifying the total numbers of rows and
% columns of the matrix, and the rest of the data are of type specified
% by [precision], which can be 'double', 'integer*4', or anyother that
% is accepted by MATLAB.  By default, the storage is assumed to be
% row-major; set [order] to 'F' to load the data in column-major.
%
% Copyright C Daniel Meliza, Z Chi 2010.  Licensed for use under Creative
% Commons Attribution-Noncommercial-Share Alike 3.0 United States
% License (http://creativecommons.org/licenses/by-nc-sa/3.0/us/).

if (nargin < 3)
  order = 'C';
end

FID=fopen(fname, 'r');
nx=fread(FID, 1, 'integer*4');
ny=fread(FID, 1, 'integer*4');

if strcmp(lower(precision),'complex')
     M=fread(FID,nx*ny*2, 'double');
     M=reshape(M,2,nx*ny)';
     M=complex(M(:,1), M(:,2));
else
     M=fread(FID, nx*ny, precision);
end
fclose(FID);
if order=='C'
  M=transpose(reshape(M, ny, nx)); 
else
  M=reshape(M,nx,ny);
end