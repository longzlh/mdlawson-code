function disp(obj)
%DISP   Display information about an RKFUNB.

  c = obj.coeffs;

  s = size(c{1},1);
  fprintf('\tRKFUNB object of block size %d-by-%d and type (%d, %d).\n', s, s, type(obj));

  if isreal(obj.K) && isreal(obj.H), f = 'Real';
  else                               f = 'Complex'; end
  s = [f '-valued Hessenberg pencil (H, K) of size '];
  fprintf('\t%s%d-by-%d.\n', s, size(obj.H, 1), size(obj.H, 2));
  
  r = 1; 
  for j = 1:length(c),
      if ~isreal(c{j}), r = 0; break; end
  end
  s = 1;
  for j = 1:length(c),
      if ~issparse(c{j}), s = 0; break; end
  end
  if r, f = 'Real';
  else  f = 'Complex'; 
  end
  if s, f = [ f ' sparse'];
  else  f = [ f ' dense']; 
  end
  
  s = [f ' coefficient matrices of size '];
  fprintf('\t%s%d-by-%d.\n', s, size(obj));
  

  
  
  if isa(c, 'sym')
    fprintf('\tVariable precision arithmetic (VPA) activated.\n');
  elseif isa(c, 'mp')
    fprintf(['\tMultiple precision arithmetic (ADVANPIX)' ...
             ' activated.\n']);
  end


end % function
