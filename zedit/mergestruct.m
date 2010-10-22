function S = mergestruct(S, X)
%
% MERGESTRUCT Merges two structures, overwriting any shared fields
% with the values stored in the second structure
%
% S = MERGESTRUCT(S, X) adds all the fields in X to S. If S has a
% field present in X, it will be replaced with the value in X.
%
% S = MERGESTRUCT(S) returns S
%
% 8/2006, CDM

if nargin > 1 && isstruct(X)
     for i = fieldnames(X)'
          fn = i{1};
          S.(fn) = X.(fn);
     end
end
