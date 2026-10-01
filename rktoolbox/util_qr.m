function [Q, R, perm] = util_qr(A, ip)
% UTIL_QR   Gram-Schmidt QR in nonstandard inner product
%
%   [Q, R, perm] = util_qr(A, ip)
% 
% computes a QR decomposition of A such that A(:,perm) = Q*R, 
% the columns of Q are orthogonal in the inner product ip(x,y),
% and the diagonal entries of R are nonincreasing.

if nargin < 2,
    ip = @(x,y) y'*x;
end

% not done: sort rows by decreasing inf-norm
% [~,ind] = sort(max(abs(A),[],2),'descend');

[m, n] = size(A); perm = 1:n;
Q = A;  
R = zeros(min(m,n),n);

% Get norms for column pivoting.
for j = 1:n,
    nrm(j) = abs(ip(Q(:,j),Q(:,j))); % Square of norms
end

% Perform 'right-affecting' MGS with column pivoting.
for j = 1:min(m,n),
    % pivoting    
    [~,p] = max(nrm(j:n)); p = j-1+p;
    Q(:,[j,p]) = Q(:,[p,j]);
    R(:,[j,p]) = R(:,[p,j]);
    perm([j,p]) = perm([p,j]);
    nrm([j,p]) = nrm([p,j]);
    
    % First orthogonalization
    R(j,j) = sqrt(ip(Q(:,j),Q(:,j)));
    if R(j,j) ~= 0,  
        Q(:,j) = Q(:,j)/R(j,j);
    end
    
    for k = j+1:n,
        R(j,k) = ip(Q(:,k),Q(:,j));
    end
    Q(:,j+1:n) = Q(:,j+1:n) - Q(:,j)*R(j,j+1:n);
    
    % Downdating norms of columns k+1:n (remaining pivot candidates)
    % abs() only because quantity may get negative due to rounding.
    nrm(j+1:n) = abs(nrm(j+1:n) - abs(R(j,j+1:n)).^2);
end

Q = Q(:,1:min(m,n));

% Reorthogonalization
R1 = zeros(min(m,n));
for j = 1:min(m,n),    
    R1(j,j) = sqrt(ip(Q(:,j),Q(:,j)));
    if R1(j,j) ~= 0,
        Q(:,j) = Q(:,j)/R1(j,j);
    end

    for k = j+1:min(m,n),
        R1(j,k) = ip(Q(:,k),Q(:,j));
    end
    Q(:,j+1:min(m,n)) = Q(:,j+1:min(m,n)) - Q(:,j)*R1(j,j+1:min(m,n));
end
R = R1*R;

