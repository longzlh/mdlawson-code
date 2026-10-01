function varargout = type(obj)
%TYPE    Return the type of a BARYFUN.

  d = length(obj.zk)-1;
  if nargout <= 1
    varargout{1} = [d, d];
  else
    varargout{1} = d;
    varargout{2} = d;
  
end