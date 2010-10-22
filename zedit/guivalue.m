function v = guivalue(obj, PN, PV)
%
% GUIVALUE Retrieves or sets a value in the guidata.
%
% If the value stored in guidata(obj) is a structure, this function
% will retrieve or set the value of the field PN in that
% structure. 
%
% V = GUIVALUE(handle, PN) Returns the value of the field PN. If
% the guidata value is not a structure, or the field does not
% exist, returns [].
%
% GUIVALUE(handle, PN, PV) Stores the value PV in the field
% PN. Raises an error if the guidata is not a structure or the
% empty matrix.
%
% 8/2006, CDM
error(nargchk(2,3,nargin))

data = guidata(obj);
v   = [];
if nargin == 2 
     % getter case
     if (~isstruct(data) || ~isfield(data,PN))
          return
     end
     v = data.(PN);
else
     % setter case
     data.(PN) = PV;
     v = PV;
     guidata(obj,data);
end
     
