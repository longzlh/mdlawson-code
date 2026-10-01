function varargout = type(obj)
%TYPE    Return the type (m+k,m) of an RKFUNB.

  k  = obj.k;
  s = size(obj.coeffs{1},1);
  mk = round(size(obj.H, 2)/s);
  
  
  if k <= 0, m = mk;
  else       m = mk - k; end
  
  if nargout <= 1
    varargout{1} = [m+k, m];
  else
    varargout{1} = m+k;
    varargout{2} = m;
  end
  
end