function [rat,pol,res,zer,z,f,w,errvec] = util_aaa(F,Z,tol,mmax,tc)
% UTIL_AAA    AAA algorithm for rational approximation.
%
% This is the AAA (adaptive Anderson-Antoulas) algorithm from [2].
% It only differs from the original algorithm in that
%
% - the output is returned in RKFUN or RKFUNM format, using the 
%   conversion described in [1],
% - it also works for a matrix-valued function F by using a scalar 
%   surrogate function f(z) = u'*F(z)*v for the AAA sampling.
%
% Usage:   [rat,pol,res,zer,z,f,w,errvec] = util_aaa(F,Z,tol,mmax,tc)
%
% Inputs:  F    = vector of scalar data values, or a function handle
%                 to a matrix-valued function
%          Z    = vector of sample points
%          tol  = relative tolerance tol, set to 1e-13 if omitted
%          mmax : type is (mmax-1,mmax-1), set to 100 if omitted
%          tc   : transform coefficients (default 1). If set to 0 the 
%                 RKFUNM coefficient transform omitted, with the consequence 
%                 that RKFUNM does no longer represent the interpolated
%                 function. However, the linearization is still valid 
%                 for nonlinear eigenvalue problems. 
%
% Outputs: rat = AAA approximant to F in RKFUN or RKFUNM format
%          pol,res,zer = vectors of poles, residues, zeros of the 
%                        scalar (surrogate) approximant
%          z,f,w = support pts, (scalar) function values, weights
%          errvec = vector of errors at each step of the AAA algorithm.
%
% WARNING: In the case of a matrix-valued function F the input error 
%          tolerance and the entries in the errvec output correspond 
%          to errors for a scalar surrogate problem f(z) = u'*F(z)*v, 
%          with random unit vectors u, v. It is up to the user to 
%          verify the approximation accuracy of the RKFUNM output  
%          for the original matrix-valued function F.
%
% References:
%
% [1] S. Elsworth and S. G{\"u}ttel, _Conversions between barycentric, 
%     RKFUN, and Newton representations of rational interpolants,_
%     Linear Algebra Appl., 576:246--257, 2019.
%
% [2] Y. Nakatsukasa, L. N. Trefethen, and O. S{\'e}te,
%     _The AAA algorithm for rational approximation,_
%     SIAM J. Sci. Comput. 40(3):A1494-A1522, 2018.

M = length(Z); % number of sample points
if nargin<3, tol = 1e-13; end % default relative tol 1e-13
if nargin<4, mmax = 100; end % default max type (99,99)
if nargin<5, tc = 1; end % get RKFUNM coefficients in b_0 = 1 basis

% treat different cases of input data
if ~isfloat(F),
    F1 = F(Z(1)); % evaluate to get dimension
    [m,n] = size(F1);
    if m == 1 && n == 1,
        f = F(Z);
        f = f(:); 
    else
        state = rng(); rng(0);
        u = randn(1,m); u = u/norm(u); v = randn(n,1); v = v/norm(v);
        rng(state);
        f = zeros(M,1);
        for j = 1:M,
            f(j) = u*F(Z(j))*v;
        end
    end
else
    f = F(:);
    m = 1; n = 1;
end
Z = Z(:);

% call AAA algorithm
[r,pol,res,zer,z,f,w,errvec] = aaa(f,Z,tol,mmax);

% convert to scalar RKFUN format
[rat,QQ] = util_bary2rkfun(z, f, w);

% in matrix-case convert to rkfunm and substitute coefficients
if m > 1 || n > 1,
    rat = rkfunm(rat);
    % first evaluate F at the sample points
    C = cell(length(z),1);
    for i = 1:length(z),
        C{i} = F(z(i));
    end
    % Apply coefficient transform (mat-vec) to rat coeffs
    % this is necessary because in RKFUNM the first basis function
    % b_0 is constant 1. The transformation may be slow, however, 
    % and is not required when solving NEPs. Hence the option to omit.
    if tc,
        rat.coeffs = util_matcell(QQ,C);
    else
        rat.coeffs = C;
    end
    
    % if m == n, append linearization pencil. we'll use the Newton 
    % interpretation of the barycentric interpolant, hence working 
    % with the original coefficients C and the Newton data.
    if m == n, 
        d = length(z)-1;
        sigma = z;
        beta = -z(2:end).*w(1:end-1)./w(2:end);
        xi = -beta.*w(2:end)./w(1:end-1);
        k = ones(d,1); h = xi;
        ii = find(xi == 0); k(ii) = 0; h(ii) = 1;
        AB = struct();
        AB.get_matrices = @() get_matrices();
        solve(); % reset cache
        AB.solve = @(nu,mu,x) solve(nu,mu,x);
        AB.multiply = @(rho,eta,x) multiply(rho,eta,x);
        AB.real = 0;
        rat.AB = AB;
    end
end

    % some private functions for the pencil construction
    function [Ad,Bd] = get_matrices()
        Ad = spdiags([sigma(1:d).*h(1:d),[0;beta(1:d-1).*h(1:d-1)]],-1:0,d,d);
        Ad = kron(Ad,speye(n));
        for jj = 1:d,
            Ad(1:n,1+n*(jj-1):n*jj) = h(d)*C{jj};
        end
        Ad(1:n, end-n+1:end) = Ad(1:n, end-n+1:end) - (h(d)*sigma(d)/beta(d))*C{end};
        Bd = spdiags([ones(d,1).*h(1:d),[0;beta(1:d-1).*k(1:d-1)]],-1:0,d,d);
        Bd = kron(Bd,speye(n));
        for jj = 1:d,
            Bd(1:n,1+n*(jj-1):n*jj) = k(d)*C{jj};
        end
        Bd(1:n, end-n+1:end) = Bd(1:n, end-n+1:end) - (h(d)/beta(d))*C{end};  
    end

    function y = multiply(rho,eta,x)
        y = zeros(size(x));
        for jj = 1:d, 
            y(1:n,:) = y(1:n,:) + (rho*h(d))*(C{jj}*x(1+(jj-1)*n:jj*n,:)) ...
                - (eta*k(d))*(C{jj}*x(1+(jj-1)*n:jj*n,:));
        end
        y(1:n,:) = y(1:n,:) - (rho*h(d)*sigma(d)/beta(d))*(C{end}*x(1+(d-1)*n:d*n,:)) ...
            + (eta*h(d)/beta(d))*(C{end}*x(1+(d-1)*n:d*n,:));
        for jj = 2:d,
            y(1+(jj-1)*n:jj*n,:) = (rho*h(jj-1)*sigma(jj-1) - eta*h(jj-1))*x(1+(jj-2)*n:(jj-1)*n,:) + ...
                (rho*h(jj-1)*beta(jj-1) - eta*k(jj-1)*beta(jj-1))*x(1+(jj-1)*n:jj*n,:);
        end
    end

    function y = solve(nu,mu,x)
        
        persistent CACHE  
        if nargin == 0, % activate/reset cache use   
            CACHE.L = {}; CACHE.U = {}; CACHE.P = {}; CACHE.Q = {}; 
            CACHE.shift = [];
            return
        end
        
        %{
        [Ad,Bd] = getMatrices();
        y_ex = (nu*Ad - mu*Bd)\x;
        % check factorization first
        M_ex = (nu*Ad - mu*Bd);
        for jj = 1:d,
            E{jj} = (nu*h(d) - mu*k(d))*C{jj};
        end
        E{d} = E{d} - (nu*h(d)*sigma(d)/beta(d) - mu*h(d)/beta(d))*C{d+1};
        vecc = nu*h.*sigma(1:d) - mu*h;
        vecd = nu*h.*beta - mu*k.*beta; 
        vecg = vecc./vecd;
        F = cell(1,d);
        F{d} = E{d};
        for jj = d:-1:2,
            F{jj-1} = E{jj-1} - vecg(jj-1)*F{jj};  
        end
        M1t = [];
        for jj = 1:d,
            M1t = [ M1t , F{jj} ];
        end
        M1 = kron(diag([0;vecd(1:d-1)]),eye(n)); 
        M1(1:n,:) = M1t;
        M2 = spdiags([vecg,ones(d,1)],-1:0,d,d);
        M2 = kron(M2,eye(n));
        norm(M1*M2 - M_ex,'fro')/norm(M_ex,'fro')        
        %}

        y = zeros(size(x));
        z = zeros(size(x));
        for jj = 2:d,
            z(1+(jj-1)*n:jj*n,:) = x(1+(jj-1)*n:jj*n,:) / ...
                (nu*h(jj-1)*beta(jj-1) - mu*k(jj-1)*beta(jj-1));
        end
        rhs = x(1:n,:);
        
        Fj = (nu*h(d) - mu*k(d))*C{d} - (nu*h(d)*sigma(d)/beta(d) - mu*h(d)/beta(d))*C{d+1};
        for jj = d:-1:2, 
            rhs = rhs - Fj*z(1+(jj-1)*n:jj*n,:);
            g = (nu*h(jj-1)*sigma(jj-1) - mu*h(jj-1)) / ...
                (nu*h(jj-1)*beta(jj-1) - mu*k(jj-1)*beta(jj-1));
            Fj = (nu*h(d) - mu*k(d))*C{jj-1} - g*Fj;
        end 
        
        %z(1:n,:) = Fj\rhs;
        
        if isfield(CACHE,'shift'), % store LU factors of Fj
            ind = find(CACHE.shift == mu/nu,1,'first');
            if isempty(ind),        % compute and store factor
                ind = length(CACHE.shift)+1;
                if issparse(Fj),
                    [CACHE.L{ind},CACHE.U{ind},CACHE.P{ind},CACHE.Q{ind}] = lu(Fj);
                else
                    [CACHE.L{ind},CACHE.U{ind},CACHE.P{ind}] = lu(Fj);
                    CACHE.Q{ind} = 1;
                end
                CACHE.shift(ind) = mu/nu;
            end
            z(1:n, :) = CACHE.Q{ind}*(CACHE.U{ind}\(CACHE.L{ind}\(CACHE.P{ind}*rhs)));
        else
            z(1:n, :) = Fj\rhs;
        end
                
        y(1:n,:) = z(1:n,:);
        for jj = 2:d,
            g = (nu*h(jj-1)*sigma(jj-1) - mu*h(jj-1)) / ...
                (nu*h(jj-1)*beta(jj-1) - mu*k(jj-1)*beta(jj-1));
            y(1+(jj-1)*n:jj*n,:) = z(1+(jj-1)*n:jj*n,:) - g*y(1+(jj-2)*n:(jj-1)*n,:);
        end
    end

end


function [r,pol,res,zer,z,f,w,errvec] = aaa(F,Z,tol,mmax)
% aaa rational approximation of data F on set Z
% [r,pol,res,zer,z,f,w,errvec] = aaa(F,Z,tol,mmax)
%
% Input: F = vector of data values, or a function handle
% Z = vector of sample points
% tol = relative tolerance tol, set to 1e-13 if omitted
% mmax: max type is (mmax-1,mmax-1), set to 100 if omitted
%
% Output: r = AAA approximant to F (function handle)
% pol,res,zer = vectors of poles, residues, zeros
% z,f,w = vectors of support pts, function values, weights
% errvec = vector of errors at each step

M = length(Z); % number of sample points
%if nargin<3, tol = 1e-13; end % default relative tol 1e-13
%if nargin<4, mmax = 100; end % default max type (99,99)
%if ~isfloat(F), F = F(Z); end % convert function handle to vector
%Z = Z(:); F = F(:); % work with column vectors
SF = spdiags(F,0,M,M); % left scaling matrix
J = 1:M; z = []; f = []; C = []; % initializations
errvec = []; R = mean(F);
for m = 1:mmax % main loop
    [~,j] = max(abs(F-R)); % select next support point
    z = [z; Z(j)]; f = [f; F(j)]; % update support points, data values
    J(J==j) = []; % update index vector
    C = [C 1./(Z-Z(j))]; % next column of Cauchy matrix
    Sf = diag(f); % right scaling matrix
    A = SF*C - C*Sf; % Loewner matrix
    [~,~,V] = svd(A(J,:),0); % SVD
    w = V(:,m); % weight vector = min sing vector
    N = C*(w.*f); D = C*w; % numerator and denominator
    R = F; R(J) = N(J)./D(J); % rational approximation
    err = norm(F-R,inf);
    errvec = [errvec; err]; % max error at sample points
    if err <= tol*norm(F,inf), break, end % stop if converged
end
r = @(zz) feval(@rhandle,zz,z,f,w); % AAA approximant as function handle
[pol,res,zer] = prz(r,z,f,w); % poles, residues, and zeros
[r,pol,res,zer,z,f,w] = ...
    cleanup(r,pol,res,zer,z,f,w,Z,F); % remove Frois. doublets (optional)
end

function [pol,res,zer] = prz(r,z,f,w) % compute poles, residues, zeros
m = length(w); B = eye(m+1); B(1,1) = 0;
E = [0 w.'; ones(m,1) diag(z)];
pol = eig(E,B); pol = pol(~isinf(pol)); % poles
dz = 1e-5*exp(2i*pi*(1:4)/4);
res = r(bsxfun(@plus,pol,dz))*dz.'/4; % residues
E = [0 (w.*f).'; ones(m,1) diag(z)];
zer = eig(E,B); zer = zer(~isinf(zer)); % zeros
end

function r = rhandle(zz,z,f,w) % evaluate r at zz
zv = zz(:); % vectorize zz if necessary
CC = 1./bsxfun(@minus,zv,z.'); % Cauchy matrix
r = (CC*(w.*f))./(CC*w); % AAA approx as vector
ii = find(isnan(r)); % find values NaN = Inf/Inf if any
for j = 1:length(ii)
    r(ii(j)) = f(find(zv(ii(j))==z)); % force interpolation there
end
r = reshape(r,size(zz)); % AAA approx
end

function [r,pol,res,zer,z,f,w] = cleanup(r,pol,res,zer,z,f,w,Z,F)
m = length(z); M = length(Z);
ii = find(abs(res)<1e-13); % find negligible residues
ni = length(ii);
if ni == 0, return, end
fprintf('%d Froissart doublets\n',ni)
for j = 1:ni
    azp = abs(z-pol(ii(j)));
    jj = find(azp == min(azp),1);
    z(jj) = []; f(jj) = []; % remove nearest support points
end
for j = 1:length(z)
    F(Z==z(j)) = []; Z(Z==z(j)) = [];
end
m = m-length(ii);
SF = spdiags(F,0,M-m,M-m);
Sf = diag(f);
C = 1./bsxfun(@minus,Z,z.');
A = SF*C - C*Sf;
[~,~,V] = svd(A,0); w = V(:,m); % solve least-squares problem again
r = @(zz) feval(@rhandle,zz,z,f,w);
[pol,res,zer] = prz(r,z,f,w); % poles, residues, and zeros
end
