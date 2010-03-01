function [X,I,J] = feat_merge(M, IND, varargin)
%
% FEAT_MERGE  Merges two or more features in a feature map
%
% X = FEAT_MERGE(M, IND)  Merges the features referred to by the
% elements of IND into a single feature. A new feature map is
% returned as X.  The new merged feature will have the same index
% as the first element in IND.  The other features maintain their
% indices.
%
% [X,I,J] = FEAT_MERGE(M, IND, 'renorm') Merges features as above, but
% renames all the features. No guarantees are made as to the
% mapping of feature numbers. The feature mapping is given by the
% returned vectors I and J, with I containing the list of the old
% feature numbers, and J the list of the new numbers.
%
% ... = FEAT_MERGE(M, IND, 'delete') causes the features to be
% deleted, rather than merged. 'delete' can be combined with 'renorm'
%
% If IND contains invalid feature numbers, they will be skipped. If
% there are no valid feature numbers in IND, or only one, X will be
% a copy of M. 
%

X = M;
I = unique(M);
J = I;

IND = IND(IND > 0 & IND <= max(I));

if ~isempty(IND)
     if ~isempty(strmatch('delete',varargin))
          replace = 0;
     else 
          replace = IND(1);
     end
     

     newfeat_ind = ismember(M,IND);
     X(newfeat_ind) = replace;
     J(ismember(I,IND)) = replace;
end


% renormalize if required
if ~isempty(strmatch('renorm',varargin))
     M = X;
     Jt = J;
     indices = unique(M(M>0));
     for i = 1:length(indices)
          X(M==indices(i)) = i;
          J(Jt==indices(i)) = i;
     end
end
