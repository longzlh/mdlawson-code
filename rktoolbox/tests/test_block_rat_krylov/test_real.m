function check = test_real()
% Tests:    - Constructing real decomposition
%           - Extending real decomposition
%           - Rerunning real decomposition

tol = 1e-13;

% Constructing real decomposition.
N = 200;
s = 5;
A = rand(N);
b = rand(N,s);
xi = [2, 1+1i, 1-1i, 4i, -4i, 0, inf];
param.real = 1;
[V,K,H] = rat_krylov(A, b, xi, param);
nrm1 = norm(A*V*K - V*H, 'fro')/norm(H, 'fro');
nrm2 = norm(V'*V - eye(size(V,2)));
r1 = [isreal(K), isreal(H)];

% Extending real decomposition
xi = [3.27 - 1.43i, 3.27 + 1.43i, 6, 1i, -1i];
param.extend = 5;
[V, K, H] = rat_krylov(A, V, K, H, xi, param);
nrm3 = norm(A*V*K - V*H, 'fro')/norm(H, 'fro');
nrm4 = norm(V'*V - eye(size(V,2)));
r2 = [isreal(K), isreal(H)];

% Rerunning real decomposition
N = 200;
s = 5;
A = rand(N, N);
b = rand(N,s);
xi = [1+3i, 1-3i, 6i, -6i, 1.3+1.42i, 1.3-1.42i ];
[V, K, H] = rat_krylov(A, b, xi, param);
nrm5 = norm(A*V*K-V*H,'fro') / norm(H, 'fro');
nrm6 = norm(V'*V-eye(size(V, 2)));

b2 = rand(N,s);
Vr = rat_krylov(A, b2, K, H);
nrm7 = norm(A*Vr*K-Vr*H,'fro')/norm(H, 'fro');

check = [[nrm1 nrm2 nrm3 nrm4 nrm5 nrm6 nrm7] < tol, r1, r2];
end