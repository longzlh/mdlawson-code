function check = test_block_polynomial()
% Tests:    - block polynomial Arnoldi, 
%           - 'CGS', 
%           - 'refinement', 
%           - matrix pencils,
%           - structures.

N = 100;
A = gallery('tridiag',N);
B = gallery('grcar', N, 3);
b = ones(N, 2);
b(1,1)   = 0;
b(4,1)   = -1;
b(N,1) = 0;
b(16,2) = -10;
b(17,2) = 3;

xi = [inf, inf, inf, inf, inf];

tol = 1e-13;

[V, K, H] = rat_krylov(A, b, xi);
nrm1 = norm(A*V*K - V*H, 'fro')/norm(H, 'fro');
nrm2 = norm(V'*V - eye(size(V,2)));

[V, K, H] = rat_krylov(A, speye(N), b, xi);
nrm3 = norm(A*V*K - V*H, 'fro')/norm(H, 'fro');
nrm4 = norm(V'*V - eye(size(V,2)));

[V, K, H] = rat_krylov(A, B, b, xi);
nrm5 = norm(A*V*K - B*V*H, 'fro')/norm(H, 'fro');
nrm6 = norm(V'*V - eye(size(V,2)));
      
param.orth = 'CGS';
param.refinement = 1;
[V, K, H] = rat_krylov(A, B, b, xi, param);
nrm7 = norm(A*V*K - B*V*H, 'fro')/norm(H, 'fro');
nrm8 = norm(V'*V - eye(size(V,2)));

[V, K, H] = rat_krylov(struct('A', A, 'B', B), b, xi);
nrm9 = norm(A*V*K - B*V*H, 'fro')/norm(H, 'fro');
nrm10 = norm(V'*V - eye(size(V,2)));

check = [nrm1 nrm2 nrm3 nrm4 nrm5 nrm6 nrm7 nrm8 nrm9 nrm10] <tol;
end 