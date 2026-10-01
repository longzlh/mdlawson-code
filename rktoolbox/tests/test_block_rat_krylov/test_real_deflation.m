function check = test_real_deflation()
% Tests:    - all columns deflate during real variant
%           - some columns deflate during real variant
%           - matching columns deflate during real variant

tol = 1e-13;

% all columns deflate during real variant
N = 200;
param.real = 1;
param.deflation_tol = 1e-10;
A  = diag(ones(N-1,1),-1);
b = zeros(N,2);
b(N-6,1) = 1;
b(N-3,2) = 1;
xi = [inf, 2+3i, 2-3i, 2];
[V, K, H, out] = rat_krylov(A, b, xi, param);
nrm1 = norm(A*V*K - V*H, 'fro')/norm(H, 'fro');
nrm2 = norm(V'*V - eye(size(V,2)));
block1 = out.blocksizes == [2, 2, 0, 0, 2];

% some columns deflate during real variant
b = zeros(N,2);
b(1,1) = 1;
b(N-2,2) = 1;
xi = [inf, 2+3i, 2-3i, 2, 7];
[V, K, H, out] = rat_krylov(A, b, xi, param);
nrm3 = norm(A*V*K - V*H, 'fro')/norm(H, 'fro');
nrm4 = norm(V'*V - eye(size(V,2)));
block2 = out.blocksizes == [2, 2, 0, 0, 2, 1];

% matching columns deflate during real variant
b = zeros(N,2);
b(1,1) = 1;
b(N-1,2) = 1;
xi = [inf, 2+3i, 2-3i, 2, 7];
[V, K, H, out] = rat_krylov(A, b, xi, param);
nrm5 = norm(A*V*K - V*H, 'fro')/norm(H, 'fro');
nrm6 = norm(V'*V - eye(size(V,2)));
block3 = out.blocksizes == [2, 2, 1, 1, 1, 1];

check = [[nrm1 nrm2 nrm3 nrm4 nrm5 nrm6] < tol, block1, block2, block3];
end