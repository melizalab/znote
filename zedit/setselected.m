function out  = setselected(handle, varargin)
%
% SETSELECTED   Sets the selected item in a list or popupmenu
%
% V = SETSELECTED(H, 'S'): Sets the selected item in the handle object H
% to S. Returns the index of the item in the String property. If
% S does not exist, no change to selection is made, and V is empty. S can
% also be a cell array of strings, in which case SETSELECTED will try to
% select all the items in S, and returns a vector of indices.
%
% S = SETSELECTED(H, V): Sets the value parameter of handle object H to V
% and returns the selected string. If V is a vector of values and the
% object is capable of multiple selections, selects all the valid indices
% in V and returns a cell array of string selections.  If no values in V
% are valid, returns [].
%
% V = SETSELECTED(H, item, choices): Sets the String property of H to
% {choices}, and selects <item>. As above, <item> can be a string, a cell
% array of strings, or a vector of indices.
%
% $Id: setselected.m,v 1.1 2006/01/31 19:59:08 meliza Exp $

if nargin < 1 || isempty(handle) || ~ishandle(handle)
    error('setselected:invalidHandle',...
        'A valid handle must be supplied as the first argument.');
elseif nargin < 2
    error('setselected:invalidArgument', ...
        'At least one choice must be supplied.');
end

if nargin > 2 && ~isempty(varargin{2})
    choices = cellwrap(varargin{2});
    set(handle, 'String', choices);
else
    choices = cellwrap(get(handle, 'String'));
end

item     = varargin{1};
maxind   = length(choices);
if isnumeric(item)
    V   = item(item > 0 & item <= maxind);
elseif ischar(item)
    V   = strmatch(item, choices, 'exact');
elseif iscell(item)
    V   = [];
    for i = 1:length(item)
        V   = cat(1,V, strmatch(item{i}, choices, 'exact'));
    end
else
    error('setselected:invalidArgument', ...
        '%s is not valid for selections of class %s', mfilename, ...
        class(item))
end
    
if ~isempty(V)
    set(handle, 'Value', V);
    if isnumeric(item)
        out     = choices{V};
    else
        out     = V;
    end
end
