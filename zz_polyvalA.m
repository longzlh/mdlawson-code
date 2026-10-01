function [y, W] = zz_polyvalA(d,H,s)
% 
% polyfitA:  Vx*c = f;  Vx = Q*R; Q*(R*c) = Q*d  =  f;
% polyvalA:  y = Vs*c = W*R*c = W*d; Vs = W*R;
% Ref: P. D. Brubeck, Y. Nakatsukasa and L. N. Trefethen,
%      Vandermonde with Arnoldi, SIAM Rev., 63 (2021), 405-415.
% 
% 
% check---------------------------------------
d = d(:); 
II = find(abs(d)>1e-16, 1,'last'); 
if isempty(II) 
    warning('Input coeffs dk are all zeros');
    y = 0;  W = []; 
    return
else
    d = d(1:II); 
    H = H(1:II,1:II-1); 
end

% -----------------------------------------



M = length(s);
W = ones(M,1) ;
n = size(H,2);
for k = 1:n
    w = s.*W(:,k);
    for j = 1:k
 
        w = w - H(j,k)*W(:,j);
 
    end
%      w = w - W(:,1:k)*H(1:k,k);
    W = [W,  w/H(k+1,k)];
end
y = W*d;
end