function check = test_extend
% Tests:    - extending, 
%           - extending deflated decomposition,
%           - extending different number of vectors.

tol = 1e-13;
deflation_tol = 1e-10;

% Extending
N = 100;
s = 6;
A = rand(N, N) + 1i*rand(N, N);
B = rand(N, N) + 1i*rand(N, N);
b = rand(N, s) + 1i*rand(N, s);
xi = [5, 7, 4+3i, 2, inf, 0, 9];
[V, K, H] = rat_krylov(A, B, b, xi);
nrm1 = norm(A*V*K - B*V*H, 'fro')/norm(H, 'fro');
nrm2 = norm(V'*V - eye(size(V,2)));

[VV, KK, HH] = rat_krylov(A, B, b, xi(1:5));
param.extend = s;
[VV, KK, HH] = rat_krylov(A, B, VV, KK, HH, xi(6:7), param);
nrm3 = norm(A*VV*KK - B*VV*HH, 'fro')/norm(H, 'fro');
nrm4 = norm(VV'*VV - eye(size(VV,2)));
nrm5 = norm(K-KK, 'fro');
nrm6 = norm(H-HH, 'fro');
nrm7 = norm(V-VV, 'fro');

% Extending deflated decomposition
param.deflation_tol = 1e-10;
param.orth = 'CGS';
param.continuation = 'last';
N = 100;
s = 4;
A = rand(N, N) + 1i*rand(N, N);
b = rand(N, s) + 1i*rand(N, s);
b = [b, b(:,1) + A^2*b(:,3), A^3*b(:,4)];
xi = rand(1,7) + 1i*rand(1,7);
[V, K, H] = rat_krylov(A, b, xi, param);
nrm8 = norm(A*V*K - V*H, 'fro')/norm(H, 'fro');
nrm9 = norm(V'*V - eye(size(V,2)));

[VV, KK, HH] = rat_krylov(A, b, xi(1:5), param);
param.extend = 6;
[VV, KK, HH] = rat_krylov(A, VV, KK, HH, xi(6:7), param);
nrm10 = norm(A*VV*KK - VV*HH, 'fro')/norm(H, 'fro');
nrm11 = norm(VV'*VV - eye(size(VV,2)));

% Extending different number of vectors
N = 100;
s = 4;
A = rand(N, N) + 1i*rand(N, N);
b = rand(N, s) + 1i*rand(N, s);
xi = rand(1,4) + 1i*rand(1,4);
[V, K, H] = rat_krylov(A, b, xi);
nrm12 = norm(A*V*K - V*H, 'fro')/norm(H, 'fro');
nrm13 = norm(V'*V - eye(size(V,2)));
[VV, KK, HH] = rat_krylov(A, V, K, H, [0, inf]);
nrm14 = norm(A*VV*KK - VV*HH, 'fro')/norm(H, 'fro');
nrm15 = norm(VV'*VV - eye(size(VV,2)));
param2.extend = 10;
[VV, KK, HH] = rat_krylov(A, V, K, H, [0, inf], param2);
nrm16 = norm(A*VV*KK - VV*HH, 'fro')/norm(H, 'fro');
nrm17 = norm(VV'*VV - eye(size(VV,2)));

check = [[nrm1 nrm2 nrm3 nrm4 nrm5 nrm6 nrm7 nrm9 nrm11 nrm12 ...
    nrm13 nrm14 nrm15 nrm16 nrm17] < tol, [nrm8 nrm10] < deflation_tol];
end