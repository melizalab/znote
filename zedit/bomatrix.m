% function bomatrix(M, fname, precision);
%
% Write a matrix into a binary file.  The first two numbers in the
% file are of integer*4 type, specifying the total numbers of rows and
% columns of the matrix, and the rest of the data are of type specified
% by [precision].  The entries are organized by row in the data.

% $Id: bomatrix.m,v 1.1 2004/05/25 00:09:01 chi Exp chi $
function bomatrix(M, fname, precision);

FID=fopen(fname, 'wb');
fwrite(FID, size(M), 'integer*4');
M=M';
if strcmp(lower(precision),'complex')
     M = [real(M); imag(M)];
     M = reshape(M,1,prod(size(M)))
end

fwrite(FID, M(:), precision);
fclose(FID);
