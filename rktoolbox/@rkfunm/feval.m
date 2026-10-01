function R = feval(obj, varargin)
%FEVAL    Evaluate RKFUNM at a scalar argument.
%
% Calling syntax: R = feval(obj, z);
%
% The rational matrix-valued function R represented by the obj 
% variable can  be evaluated pointwise at a scalar argument z. 
%
% This function uses the QR evaluation of an RKFUN described on p. 18 of
%
% [1] M. Berljafa and S. G{\"u}ttel. Generalized rational Krylov
%     decompositions with an application to rational approximation,
%     SIAM J. Matrix Anal. Appl., 36(2):894--916, 2015.

z = varargin{1};
if ~isnumeric(z) || norm(size(z) - [1,1],inf)>0,
    error('FEVAL: RKFUNM evaluation only possible for scalars arguments.');
end
if length(varargin) >= 2,
   Nmax = varargin{2};
else
   Nmax = length(obj.coeffs)-1;
end

v = null((z*obj.K - obj.H)')'; 
v = v/v(1);
R = obj.coeffs{1};
for j = 2:Nmax+1, 
     R = R + v(j)*obj.coeffs{j};
end

end

