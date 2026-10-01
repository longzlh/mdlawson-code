function check = test_deflation()
% Tests:    - Arnoldi deflation, 
%           - invariance,
%           - rational deflation.

tol = 1e-12;
deflation_tol = 1e-9;

% Arnoldi deflation
param.deflation_tol = deflation_tol;
N = 200;
A = rand(N, N) + 1i*rand(N,N);
b = rand(N, 5) + 1i*rand(N, 5);
b = [b(:,2)+6*A*b(:,4), b, A^3*b(:,1)+ b(:,5), A*b(:,3)];
xi = rand(1,5) + 1i*rand(1,5);
[V, K, H, out] = rat_krylov(A, b, xi, param);
nrm1 = norm(A*V*K - V*H, 'fro')/norm(H, 'fro');
nrm2 = norm(V'*V - eye(size(V,2)));
block1 = [8,6,6,5,5,5] == out.blocksizes;

N = 200;
param.deflation_tol = deflation_tol;
A = sparse(N,N);
A(2:N+1:end) = 1;
A(1,N) = 0;
b = zeros(N,3);
b(1,1) = 1;
b(5,2) = 1;
b(N-2,3) = 1; 
xi = [1, 1, 1, 1, 1];

param.continuation = 'last';
[V, K, H, out] = rat_krylov(A, b, xi, param);
nrm3 = norm(A*V*K - V*H, 'fro')/norm(H, 'fro');
nrm4 = norm(V'*V -eye(size(V, 2)));
block2 = (out.blocksizes == [3,3,3,2,1,1]);

% Rational deflation
param.deflation_tol = deflation_tol;
N = 200;
s = 5;
A = rand(N, N) + 1i*rand(N,N);
b = rand(N, s) + 1i*rand(N, s);
xi = rand(1,4) + 1i*rand(1,4);
m = size(xi, 2);

param.continuation = 'last';
[~, K, H] = rat_krylov(A, b, xi, param);
r = eig(H(1:m*s, 1:m*s)/K(1:m*s, 1:m*s));
xi = [xi, r(1)];
[V, K, H] = rat_krylov(A, b, xi, param);
nrm5 = norm(A*V*K - V*H, 'fro')/norm(H, 'fro');
nrm6 = norm(V'*V - eye(size(V,2)));

% Invariance reached
N = 200;
A = rand(N, N) + 1i*rand(N,N);
b = rand(N, 110) + 1i*rand(N, 110);
xi = rand(1,2) + 1i*rand(1,2);
[V, K, H] = rat_krylov(A, b, xi);
nrm7 = norm(A*V*K - V*H, 'fro')/norm(H, 'fro');
nrm8 = norm(V'*V - eye(size(V,2)));

check = [[nrm1 nrm2 nrm3 nrm4 nrm5 nrm6] < deflation_tol, [nrm7 nrm8] < tol, block1, block2];
end
