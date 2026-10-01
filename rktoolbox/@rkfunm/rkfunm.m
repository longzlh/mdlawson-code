%RKFUNM    RKFUNM constructor.
%
% Usages:
% 1) obj = rkfunm(K, H, coeffs, k) constructs an RKFUNM from the
%    recursion matrices (H,K), and the coefficient matrices 
%    provided in coeffs.
%    The optional parameter k specifies sub- or superdiagonal
%    rational functions.
% 2) obj = rkfun(str, params) returns an RKFUNM from the gallery. 
%    To see a list of the gallery functions type help RKFUNM/gallery.
%    Example: r = rkfun('step'); ezplot(r, [-2, 2])


classdef rkfunm
properties
K
H
coeffs
k
AB
end % properties
    
methods(Access = public, Static = true)
% Constructor. (in methods)
function obj = rkfunm(varargin)
  
   % Convert from RKFUN.
  if isa(varargin{1},'rkfun')
    obj.K = varargin{1}.K;
    obj.H = varargin{1}.H;
    obj.coeffs = num2cell(varargin{1}.coeffs);
    obj.k = varargin{1}.k;
    return
  end
    
  % Gallery.
  if ischar(varargin{1})
    obj = rkfunm.gallery(varargin{:});
    return
  end
  
  % Construct rkfun from pencil and coeffs.
  obj.K = varargin{1};
  obj.H = varargin{2};
  obj.coeffs = varargin{3};
  if nargin >= 4
    obj.k = varargin{4};
  else
    obj.k = 0;
  end
  if nargin >= 5, % hold pencil structure
    obj.AB = varargin{5};
  else
    obj.AB = [];
  end
  
end % function rkfunm
        
function obj = gallery(varargin)
% GALLERY    Collection of rational matrix-valued functions.
%
% obj = rkfunm.gallery(funname, param1, param2, ...) takes 
% funname, a case-insensitive string that is the name of 
% a rational function family, and the family's input 
% parameters. 
%
% At the moment the gallery is empty. 

  obj = gallery(varargin{:});

end

end % methods    
end % classdef rkfun
