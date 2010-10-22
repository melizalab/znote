function ch = kids(handle, param, val)
%
% KIDS Returns the handles of the children of a graphical handle
% object.
%
% H = KIDS(HANDLE) returns the results of get(handle,'children')
%
% H = KIDS(HANDLE, PARAM, VAL) returns the previous handles passed
% through findobj(H,'flat',PARAM,VAL)
%
%
ch = get(handle,'children');
if nargin > 1
     ch = findobj(ch, 'flat', param, val);
end
