function check = test_varargin()
% Test:     - all possible varargin

tol = 1e-10;

% all possible varargin
N = 100;
s = 5;
A = rand(N,N) + 1i*rand(N,N);
B = rand(N,N) + 1i*rand(N,N);
b = rand(N,s) + 1i*rand(N,s);
xi = rand(1,5) + 1i*rand(1,5);
AB = struct('A', A, 'B', B);

[V, K, H] = rat_krylov(A, b, xi);
nrm1 = norm(A*V*K - V*H, 'fro')/norm(H, 'fro');
nrm2 = norm(V'*V - eye(size(V,2)));
[VV, KK, HH] = rat_krylov(A, V, K, H, xi);
nrm3 = norm(A*VV*KK - VV*HH, 'fro')/norm(HH, 'fro');
nrm4 = norm(VV'*VV - eye(size(VV,2)));
[V, K, H] = rat_krylov(A, b, K, H);
nrm5 = norm(A*V*K - V*H, 'fro')/norm(H, 'fro');

[V, K, H] = rat_krylov(A, B, b, xi); 
nrm6 = norm(A*V*K - B*V*H, 'fro')/norm(H, 'fro');
nrm7 = norm(V'*V - eye(size(V,2)));
[VV, KK, HH] = rat_krylov(A, B, V, K, H, xi);
nrm8 = norm(A*VV*KK - B*VV*HH, 'fro')/norm(HH, 'fro');
nrm9 = norm(VV'*VV - eye(size(VV,2)));
[V, K, H] = rat_krylov(A, B, b, K, H);
nrm10 = norm(A*V*K - B*V*H, 'fro')/norm(H, 'fro');


[V, K, H] = rat_krylov(AB, b, xi) ;
nrm11 = norm(A*V*K - B*V*H, 'fro')/norm(H, 'fro');
nrm12 = norm(V'*V - eye(size(V,2)));
[VV, KK, HH] = rat_krylov(AB, V, K, H, xi);
nrm13 = norm(A*VV*KK - B*VV*HH, 'fro')/norm(HH, 'fro');
nrm14 = norm(VV'*VV - eye(size(VV,2)));
[V, K, H] = rat_krylov(AB, b, K, H);
nrm15 = norm(A*V*K - B*V*H, 'fro')/norm(H, 'fro');


check = [nrm1 nrm2 nrm3 nrm4 nrm5 nrm6 nrm7 nrm8 nrm9 ...
    nrm10 nrm11 nrm12 nrm13 nrm14 nrm15] < tol;
end