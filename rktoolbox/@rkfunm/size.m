function s = size(obj, dim)
%SIZE    Return the size of an RKFUNM.
  
  s = size(obj.coeffs{1});
  if nargin == 2,  s = s(dim); end

end