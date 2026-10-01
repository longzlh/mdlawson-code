function check = test_qr_permutation()
% Test:     - rank revealing qr
%           - column premutations during deflation

tol = 1e-12;

% Rank revealing QR;
I = eye(6);
ip = @(x, y) y'*x;
A = [I(:,1), I(:,2), I(:,1) + I(:,2), I(:,3), I(:,4)];

[Q, R, E] = util_qr(A, ip); 
Einv = E;
Einv(E) = 1:length(E);
nrm1 = norm(A(:,E) - Q*R);

permQ = Q(:, Einv);
permR = R(Einv, Einv);
nrm2 = norm(A - permQ*permR);

r = 4;
newQ = permQ;
newQ(:, E(r+1:5)) = [];
newR = permR;
newR(E(r+1:5),:) = [];

nrm3 = norm(A - newQ*newR);


% Column permutation during deflation
N = 50;
A = gallery('tridiag', N);

I = eye(N);
b = [I(:,1), I(:,3), A*I(:,3), I(:,10), I(:,11)];

xi = [inf];
param.deflation_tol = 1e-10;

[V, K, H, out] = rat_krylov(A, b, xi, param);
nrm4 = norm(A*V*K - V*H, 'fro')/norm(H, 'fro');
nrm5 = norm(V'*V - eye(size(V,2)));

deflated_columns = out.column_deflation;

check = [[nrm1 nrm2 nrm3 nrm4 nrm5] < tol, deflated_columns == [1, 1, 0, 1, 1]];
end