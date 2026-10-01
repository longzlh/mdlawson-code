function xi = poles(obj)
%POLES    Return the poles of an RKFUNB.
  
  s = size(obj.coeffs{1},1); % block size
  xi = util_pencil_poles(obj.K, obj.H, s);
  xi = xi(:);
  
  % Only output the m smallest poles, where m is denominator
  % degree.
  t = type(obj);
  [~, ind] = sort(abs(xi), 'ascend');
  xi = xi(ind);
  xi = xi(1:t(2));
  
end
