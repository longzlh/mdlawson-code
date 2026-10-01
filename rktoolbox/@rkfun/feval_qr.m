function F = feval_qr(r,Z)
%FEVAL_QR    Point-wise evaluation of an RKFUN.
% This is an alternative to evaluation via rerunning and can
% achieve better results when the K-matrix in the RKFUN is
% badly conditioned [1].
%
% Calling syntax: w = feval_qr(obj, z);
%
% [1] M. Berljafa and S. G{\"u}ttel. Generalized rational Krylov
%     decompositions with an application to rational approximation,
%     SIAM J. Matrix Anal. Appl., 36(2):894--916, 2015.

[m,n] = size(Z);
Z = Z(:); F = zeros(length(Z),1);
for j = 1:length(Z),
       v = null((Z(j)*r.K - r.H).').'; 
       v = v/v(1);
       F(j) = v*r.coeffs;
end
F = reshape(F,m,n);
