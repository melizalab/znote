function out = getselected(handle)
% GETSELECTED Returns the string selected in a list or popupmenu
% uicontrol. 
%
% If multiple items are selected, a character array is returned. If the
% handle refers to a uibuttongroup object, then the tag of the currently
% selected button is returned.
%
% selected = GETSELECTED(handle) 
%
% handle - graphics object handle for the list or uicontrol
%
% $Id: GetSelected.m,v 1.5 2006/01/20 00:04:45 meliza Exp $

type    = get(handle,'type');
switch type
    case 'uicontrol'
        i = get(handle,'Value');
        s = get(handle,'String');
        if any(i == 0) || isempty(s)
            out = '';
        else
            if iscell(s)
                out = char({s{i}});
            elseif i <= length(s)
                out = s(i,:);
            else
                out = s;
            end
        end
        if strcmpi(out, ' ')
            out = '';
        end
    case 'uipanel'
        % uibuttongroup is actually type 'uipanel', so we have to check if
        % the SelectedObject property exists
        if isprop(handle,'SelectedObject')
            selected    = get(handle, 'SelectedObject');
            if ishandle(selected)
                out     = get(selected, 'tag');
            else
                out     = [];
            end
        else
            error('invalidObjectType',...
                'GETSELECTED can only be used with uicontrols and uibuttongroups');
        end
end
            
