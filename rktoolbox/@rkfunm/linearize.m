function AB = linearize(varargin)
%LINEARIZE    Linearizes an RKFUNM and returns a pencil structure 
%             to be used by RAT_KRYLOV.
% This is the same as UTIL_LINEARISE_NLEP (with an S)!

% Call the English function (has better manners):
  AB = linearise(varargin{:});
  
end