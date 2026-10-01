function [rat,QQ] = util_bary2rkfun(z, f, w)
% UTIL_BARY2RKFUN    Constructs an RKFUN from a barycentric interpolant.
% 
% rat = util_bary2rkfun(z, f, w) takes the vectors z, f, w all of the 
% same length m and contructs an RKFUN rat which interpolates the 
% values f at the points z. More precisely, rat satisfies
%
% rat(zz) = sum( f.*w./(zz - z) ) / sum( w./(zz - z) )
%
% This rational interpolant will be of type [m-1,m-1].
% 
% The conversion is described in 
%
% S. Elsworth and S. G{\"u}ttel, _Conversions between barycentric, RKFUN, 
% and Newton representations of rational interpolants,_ 
% MIMS Eprint 2017.xx, (<http://eprints.maths.manchester.ac.uk/2593/>), 
% Manchester Institute for Mathematical Sciences, 
% The University of Manchester, UK, 2017.

z = z(:); f = f(:); w = w(:); 
m = length(z); 
H = full(spdiags([ w(1:m-1).*z(2:m) , -w(2:m).*z(1:m-1) ],-1:0,m,m-1));
K = full(spdiags([ w(1:m-1) , -w(2:m) ],-1:0,m,m-1));

coeffs = f;

% bring const = 1 to first column of V (1 = sum of Lagrange basis funs)
[Q1,~] = qr(ones(m,1)); 
s = Q1(1,1);    
% V --> V*(Q/s)
H1 = s*(Q1'*H);
K1 = s*(Q1'*K);
coeffs1 = s*(Q1'*coeffs);
if isreal(H1) && isreal(K1), 
    [HH, KK, Q2, Z] = qz(H1(2:end,:),K1(2:end,:),'real');
else
    [HH, KK, Q2, Z] = qz(H1(2:end,:),K1(2:end,:));
end
%H2 = blkdiag(1,Q)*H1*Z;
%K2 = blkdiag(1,Q)*K1*Z;
H2 = [ H1(1,:)*Z ; HH ]; 
K2 = [ K1(1,:)*Z ; KK ]; 
coeffs2 = blkdiag(1,Q2)*coeffs1;
rat = rkfun(K2, H2, coeffs2);
QQ = blkdiag(1,Q2)*(s*Q1'); % cumulated transform of coefficients

