function check = test_inner_product()
% Tests:    - alternative inner product, 
%           - alternative bilinear form, 
%           - vector ordering during QR.

tol = 1e-12;

% alterative inner product
N = 200;
A = rand(N, N) + 1i*rand(N,N);
b = rand(N, 5) + 1i*rand(N, 5);
xi = rand(1,5) + 1i*rand(1,5);

D = diag(rand(N,1));
param.inner_product = @(x, y) y'*D*x;
[V, K, H] = rat_krylov(A, b, xi, param);
nrm1 = norm(A*V*K - V*H, 'fro')/norm(H, 'fro');
nrm2 = norm(param.inner_product(V,V) - eye(size(V,2)));

% alternative bi-linear form
N = 200;
A = rand(N, N) + 1i*rand(N,N);
b = rand(N, 5) + 1i*rand(N, 5);
xi = rand(1,5) + 1i*rand(1,5);

D = diag([rand(180,1); zeros(20,1)]);
param.inner_product = @(x, y) y'*D*x;
[V, K, H] = rat_krylov(A, b, xi, param);
nrm3 = norm(A*V*K - V*H, 'fro')/norm(H, 'fro');
nrm4 = norm(param.inner_product(V,V) - eye(size(V,2)));

% QR test: usually, the 2 vector deflates
N = 30;
b = rand(N,3);
b = [b(:,1), b(:,1)+b(:,3), b(:,2), b(:,3)];
inner_product = @(x, y) y'*x;
[w, RR, E] = util_qr(b, inner_product);
w = w(:,E);
Einv = E;
Einv(E) = 1:length(E);
RR = RR(E,Einv);
r = 3;
RR(Einv(r+1:size(b,2)),:) = [];
w(:,Einv(r+1:size(b,2))) = [];
ind = sort(E(size(b,2)-r+1:size(b,2)));
nrm5 = norm(b - w*RR);
inds = (ind == [1,3,4]);

% QR test: unexpectedly, not always the fourth vector deflates, 
N = 30;
b = rand(N,3);
b = [b, zeros(N,1)];
inner_product = @(x, y) y'*x;
[w, RR, E] = util_qr(b, inner_product);
w = w(:,E);
Einv = E;
Einv(E) = 1:length(E);
RR = RR(E,Einv);
r = 3;
RR(Einv(r+1:size(b,2)),:) = [];
w(:,Einv(r+1:size(b,2))) = [];
ind = sort(E(size(b,2)-r+1:size(b,2)));
nrm6 = norm(b - w*RR);

check = [[nrm1 nrm2 nrm3 nrm4 nrm5 nrm6] < tol, inds];
end