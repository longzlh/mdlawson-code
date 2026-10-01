function check = test_rerun
% Test:     - rerunning undeflated decomposition
%           - rerunning decomposition with pencil (A,B)

tol = 1e-12;

% rerunning undeflated decomposition.
N = 200;
s = 5;
A = rand(N, N) + 1i*rand(N, N);
b = rand(N,s) + 1i*rand(N, s);
xi = 1*rand(1, 7) + 1i*rand(1,7);

[V, K, H] = rat_krylov(A, b, xi);
nrm1 = norm(A*V*K-V*H,'fro') / norm(H, 'fro');
nrm2 = norm(V'*V-eye(size(V, 2)));

[w,~,E] = util_qr(b, @(x,y) y'*x);
bb = w(:,E);
Vr = rat_krylov(A, bb, K, H);
nrm3 = norm(V*(V'*Vr)-Vr, 'fro')/norm(Vr, 'fro');
nrm4 = norm(A*Vr*K-Vr*H,'fro')/norm(Vr, 'fro');

b2 = rand(N,s) + 1i*rand(N,s);
Vrr = rat_krylov(A, b2, K, H);
nrm5 = norm(A*Vrr*K-Vrr*H,'fro')/norm(Vrr, 'fro');


% rerunning undeflated decomposition with pencil (A, B).
N = 200;
s = 5;
A = rand(N, N);
B = rand(N, N);
b = rand(N,s);
xi = rand(1, 7) + 1i*rand(1,7);

[V, K, H] = rat_krylov(A, B, b, xi);
nrm6 = norm(A*V*K-B*V*H,'fro') / norm(H, 'fro');
nrm7 = norm(V'*V-eye(size(V, 2)));

[w,~,E] = util_qr(b, @(x,y) y'*x);
Einv = E; Einv(E) = 1:length(E);
bb = w(:,Einv);

Vr = rat_krylov(A, B, bb, K, H);
nrm8 = norm(V*(V'*Vr)-Vr)/length(xi);
nrm9 = norm(A*Vr*K-B*Vr*H,'fro')/norm(H, 'fro');

b2 = rand(N,s) + 1i*rand(N,s);
Vrr = rat_krylov(A, B, b2, K, H);
nrm10 = norm(A*Vrr*K-B*Vrr*H,'fro')/norm(Vrr, 'fro');


check = [nrm1 nrm2 nrm3 nrm4 nrm5 nrm6 nrm7 nrm8 nrm9 nrm10] < tol;
end