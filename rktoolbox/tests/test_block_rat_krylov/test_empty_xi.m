function check = test_empty_xi()
% Tests:    - construction of space with no poles, 
%           - rerunning the space, 
%           - extending the space with and without poles.  

tol = 1e-13;

% Construction
A = rand(100,100) + 1i*rand(100,100);
b = rand(100,10) + 1i*rand(100,10);
xi = [];
[V, K, H] = rat_krylov(A, b, xi);
nrm1 = norm(A*V*K - V*H, 'fro');
nrm2 = norm(V'*V - eye(size(V, 2)));

% Rerunning
b = rand(100,10) + 1i*rand(100,10);
[V, K, H] = rat_krylov(A, b, K, H);
nrm3 = norm(A*V*K - V*H, 'fro');

% Extending xi = []
xi = [];
param.extend = 10;
[VV, KK, HH] = rat_krylov(A, V, K, H, xi, param);
nrm4 = norm(A*VV*KK - VV*HH, 'fro');
nrm5 = norm(VV'*VV - eye(size(VV, 2)));
nrm6 = norm(H-HH, 'fro');
nrm7 = norm(K-KK, 'fro');

% Extending xi = [1+2i, 6, -7i]
xi = [1+2i, 6, -7i];
[V, K, H] = rat_krylov(A, V, K, H, xi);
nrm8 = norm(A*V*K - V*H, 'fro')/norm(H, 'fro');
nrm9 = norm(V'*V - eye(size(V, 2)));

check = [nrm1 nrm2 nrm3 nrm4 nrm5 nrm6 nrm7 nrm8 nrm9] < tol;
end
