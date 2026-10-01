function [V, K, H, param] = rat_krylov(varargin)
% RAT_KRYLOV    (Block) Rational Krylov method.
%
% This function is at the core of the Rational Krylov Toolbox,
% released as part of [1]. There are three main functionalities:
%   (I)   run the (block) rational Arnoldi algorithm [4, 5, 6];
%   (II)  extend a (block) rational Arnoldi decomposition (RAD);
%   (III) rerun recursion given by an RAD [2].
%
% Each of the functionalities can be invoked for:
%   (i)   a single matrix A;
%   (ii)  a matrix pencil (A, B);
%   (iii) a matrix pencil (A, B) represented by a structure AB.
%
% Additionally, a structure param may be provided to specify
% options detailed below. Instead of param, the string 'real' may be
% provided, which is equivalent to providing the structure param
% with a single field 'real' being set to 1.
%
% We now list the possible calls. The output is always of the form
% [V, K, H, param] and is therefore excluded from the following list.
%
%   rat_krylov(A, b, xi)                 % (I) - (i)
%   rat_krylov(A, b, xi, param)          % (I) - (i)
%   rat_krylov(A, V, K, H, xi)           % (II) - (i)
%   rat_krylov(A, V, K, H, xi, param)    % (II) - (i)
%   rat_krylov(A, V, K, H)               % (III) - (i)
%   rat_krylov(A, V, K, H, param)        % (III) - (i)
%
%   rat_krylov(A, B, b, xi)              % (I) - (ii)
%   rat_krylov(A, B, b, xi, param)       % (I) - (ii)
%   rat_krylov(A, B, V, K, H, xi)        % (II) - (ii)
%   rat_krylov(A, B, V, K, H, xi, param) % (II) - (ii)
%   rat_krylov(A, B, V, K, H)            % (III) - (ii)
%   rat_krylov(A, B, V, K, H, param)     % (III) - (ii)
%
%   rat_krylov(AB, b, xi)                % (I) - (iii)
%   rat_krylov(AB, b, xi, param)         % (I) - (iii)
%   rat_krylov(AB, V, K, H, xi)          % (II) - (iii)
%   rat_krylov(AB, V, K, H, xi, param)   % (II) - (iii)
%   rat_krylov(AB, V, K, H)              % (III) - (iii)
%   rat_krylov(AB, V, K, H, param)       % (III) - (iii)
%
% The input arguments are:
%
%   - A, B,   N-by-N matrices;
%   - AB,     struct representing an N-by-N pencil, see below;
%   - b,      N-by-s matrix (s column vectors);
%   - V,      N-by-ns matrix;
%   - xi,     1-by-m row-vector or empty array;
%   - K, H,   upper-Hessenberg matrices with one block-row more than column;
%   - param,  parameter structure (see below).
%
% The fields of the param structure are listed below. The names are
% self-explanatory, and if not provided, default values are used.
%
% param.deflation_tol       eps(1) (default), tolerence for detecting
%                           rank deficiency in the (block) Krylov space;
% param.orth,               'MGS' (default) or 'CGS' orthogonalization;
% param.reorth,             0, 1 (default);
% param.inner_product,      function handle of the form
%                           @(x, y) y'*x; (default);
% param.real,               0 (default), 1, 2;
%                           whether to use the real version [3] of
%                           rational Arnoldi or not;
%                           option 2 is used only for (III) with
%                           complex data on a quasi-RAD;
% param.refinement,         0 (default), 1;
%                           whether to do one step of iterative
%                           refinement for the linear system solves or
%                           not;
% param.waitbar,            0 (default), 1;
%                           whether to include a waitbar showing the
%                           progress of RAT_KRYLOV or not;
% param.continuation,       'ruhe' (default),
%                           'almost-last',
%                           'last',
%                           'near-optimal'; specifies continuation strategy
%                           to be used (currently only in non-block case);
% param.continuation_root,  can be used to specify the continuation root, see
%                           [3, Sec. 2] for building the rational Krylov
%                           space. If not set, used defaults as listed in [2].
% param.continuation_solve, function handle of the form (default)
%                           @(AB, nu, mu, x, param) AB.solve(nu, mu, x);
% param.continuation_m,     4 (default); integer which may be used
%                           by param.continuation_solve;
% param.continuation_bounds,0 (default), 1; whether to compute the norms of
%                           \widehat\vf_{j+1} and \widehat\ve_{j+1} used to
%                           produce the bounds related to [3, eq. (3.18)]
%                           and [3, eq. (3.20)]. Works only if param.p = 1,
%                           and param.continuation = 'near-optimal';
% param.extend              1 (default), the number of vectors used when
%                           extending a decomposition.
% param.balance             1 (default), 0;
%                           whether to apply column scaling to pencil or
%                           not;
%
% The following param field is used to simulate the parallel execution with p
% parallel processors, in the sense that there is less freedom for choosing the
% continuation vectors; cf. [3, Figure 4.2]. It can thus be used to assist
% further research, for instance, when designing approximate linear system
% solvers for the near-optimal continuation strategy for a particular
% application.
%
% param.p,                   1 (default); integer used to simulate number
%                            of parallel processors.
%                            Currently only available when s = 1.
%
% If a matrix pencil (A, B) is provided by the structure AB,
% type (iii) for calling the software, it should provide the following.
%
% AB.multiply,  function handle which takes as arguments
%               @(rho, eta, x) and returns rho*A*x-eta*B*x;
% AB.solve,     function handle which takes as arguments
%               @(nu, mu, x) and returns (nu*A-mu*B)\x;
% AB.isreal,    true/false flag specifying whether the pencil is
%               real-valued or not.
%
% Output arguments are:
%
%   - V,      N-by-(m+n)s matrix spanning the rational Krylov space;
%   - K, H,   (m+n)s-by-(m+n-1)s upper-Hessenberg matrices.
%   - param,  the structure from above with some additional
%             information, of interest to the developers.
%
% References:
%
% [1] M. Berljafa and S. G{\"u}ttel. Generalized rational Krylov
%     decompositions with an application to rational approximation,
%     SIAM J. Matrix Anal. Appl., 36(2):894--916, 2015.
%
% [2] M. Berljafa and S. G{\"u}ttel. The RKFIT algorithm for
%     nonlinear rational approximation, SIAM J. Sci. Comput.,
%     39(5):A2049--A2071, 2017.
%
% [3] M. Berljafa and S. G{\"u}ttel. Parallelization of the rational
%     Arnoldi algorithm, SIAM J. Sci. Comput., 39(5):S197--S221, 2017.
%
% [4] S. Elsworth and S. Guettel.
%     _The block rational Arnoldi method,_
%     MIMS Eprint 2019.2 (<http://eprints.maths.manchester.ac.uk/2685/>),
%     Manchester Institute for Mathematical Sciences, 
%     The University of Manchester, UK, 2019.
%
% [5] A. Ruhe. Rational Krylov: A practical algorithm for large sparse
%     nonsymmetric matrix pencils, SIAM J. Sci. Comput., 19(5):1535--1551,
%     1998.
%
% [6] A. Ruhe. The rational Krylov algorithm for nonsymmetric eigenvalue
%     problems. III: Complex shifts for real matrices, BIT,
%     34:165--176, 1994.

warning off backtrace
[AB, V, K, H, N, j, m, s, rerun, param] = parse_argin(varargin{:});

if param.continuation_bounds && param.p == 1
    W    = zeros(N, m+1);
    Fhat = zeros(1, m);
    fhat = zeros(1, m);
end

R = [zeros(size(K,1),s) 0*K];
T = 0*K;
xi  = param.xi;
nu  = param.moebius(1, :);
mu  = param.moebius(2, :);
rho = param.moebius(3, :);
eta = param.moebius(4, :);

realopt = param.real;

if rerun
    rerun_krylov;
else
    try
        run_krylov; 
    catch e
        error(e.message())
    end
end

param.R = R;
param.T = T;

% remove row and column indexing
param = rmfield(param,'r_b');
param = rmfield(param,'c_b');

param.moebius(1, :) = nu;
param.moebius(2, :) = mu;
param.moebius(3, :) = rho;
param.moebius(4, :) = eta;

if param.continuation_bounds && param.p == 1
    W(:, 1)    = V(:, 1);
    param.W    = W;
    param.Fhat = Fhat;
    param.fhat = fhat;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Nested functions. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function run_krylov
        bd1 = false; % Deflation occured, some (possibly all) columns deflate.
        bd2 = false; % Deflation during real and complex part of real variant does not match.
        
        offset = j-1; % If expanding the array of poles is short.
        
        % Starting vector.
        if j == 1 % If creating a new RAD.
            r_b = zeros(1, m+2); c_b = zeros(1,m+1); % Create block indexing.
            r_b(1) = 1; c_b(1) = 1;
            % Since starting block vector full rank, no deflation test needed.
            [w,RR,E] = util_qr(V(:,1:s), param.inner_product);
            Einv = E; Einv(E) = 1:length(E);
            V(:,1:s) = w(:,Einv);
            R(1:s,1:s) = RR(Einv,Einv);
            r_b(2) = s+1;
        else
            r_b = param.r_b; c_b = param.c_b;
        end
        
        if param.waitbar, hwb = waitbar(0, 'rat\_krylov'); end
        while j <= m
            
            % Compute continuation vector.
            continuation_pair;
            w = V(:, 1:r_b(j+1)-1)*T(1:r_b(j+1)-1, c_b(j):c_b(j)+s-1);
            % Polynomial part.
            w = AB.multiply(rho(j-offset), eta(j-offset), w);
            
            % Rational part, with eventual refinement.
            if param.refinement
                ww = w;
            end
            w = AB.solve(nu(j-offset), mu(j-offset), w);
            
            if param.refinement
                ww = ww - AB.multiply(nu(j-offset), mu(j-offset), w);
                w = w + AB.solve(nu(j-offset), mu(j-offset), ww);
            end
            
            if param.continuation_bounds
                W(:, j+1) = w;
            end
            
            % Orthogonalization.
            if isreal(xi(j-offset)) || realopt == 0
                switch param.orth
                    case 'CGS'
                        hh = param.inner_product(w, V(:, 1:r_b(j+1)-1));
                        w = w - V(:, 1:r_b(j+1)-1)*hh;
                        H(1:r_b(j+1)-1, c_b(j):c_b(j)+s-1) = hh;
                        if param.reorth == 1
                            hh = param.inner_product(w, V(:, 1:r_b(j+1)-1));
                            w = w - V(:, 1:r_b(j+1)-1)*hh;
                            H(1:r_b(j+1)-1, c_b(j):c_b(j)+s-1) = H(1:r_b(j+1)-1, c_b(j):c_b(j)+s-1) + hh;
                        end
                        
                        [w, RR, E] = util_qr(w, param.inner_product); % Rank revealing qr.
                        Einv = E; Einv(E) = 1:length(E);
                        w = w(:,Einv);
                        RR = RR(Einv,Einv);
                        % 'CGS' end
                        
                    case 'MGS'
                        for col = 1:j
                            hh = param.inner_product(w, V(:, r_b(col):r_b(col+1)-1));
                            w = w - V(:, r_b(col):r_b(col+1)-1)*hh;
                            H(r_b(col):r_b(col+1)-1, c_b(j):c_b(j)+s-1) = H(r_b(col):r_b(col+1)-1, c_b(j):c_b(j)+s-1) + hh;
                        end
                        if param.reorth == 1
                            for col = 1:j
                                hh = param.inner_product(w, V(:, r_b(col):r_b(col+1)-1));
                                w = w - V(:, r_b(col):r_b(col+1)-1)*hh;
                                H(r_b(col):r_b(col+1)-1, c_b(j):c_b(j)+s-1) = H(r_b(col):r_b(col+1)-1, c_b(j):c_b(j)+s-1) + hh;
                            end
                        end
                        [w, RR, E] = util_qr(w, param.inner_product); % Rank revealing qr.
                        Einv = E; Einv(E) = 1:length(E);
                        w = w(:,Einv);
                        RR = RR(Einv,Einv);

                        % 'MGS' end
                        
                end % switch param.orth
                
                % inf/nan check
                ssum = sum(sum([H(:,c_b(j):c_b(j)+s-1);RR]));
                if isinf(ssum) || isnan(ssum)
                    ME = MException('VerifyOutput:OutOfBounds', ...
                        ['Error: inf/nan occurred during orthogonalisation. Maybe pole ', ...
                       num2str(xi(j)), ' too close to eigenvalue of pencil?']);
                    throw(ME);
                end
                
                % Deflation test
                diagr = svd(RR);
                r = find(diagr > max(diagr(1)*param.deflation_tol, ...
                    norm(H(1:r_b(j+1)-1, c_b(j):c_b(j)+s-1),'fro')*param.deflation_tol), 1, 'last');
                if isempty(r), r = 0; end
                if r ~= s, bd1 = true; end
                
                r_b(j+2) = r_b(j+1) + r; c_b(j+1) = c_b(j) + s;
                RR(E(r+1:s),:) = [];
                H(r_b(j+1):r_b(j+2)-1,c_b(j):c_b(j)+s-1) = RR;
                w(:,E(r+1:s)) = [];
                V(:,r_b(j+1):r_b(j+2)-1) = w;
                inds = sort(c_b(j)-1 + E(1:r));
                param.column_deflation(inds) = ones(1,length(inds));
                
                % Setting the decomposition.
                R(1:r_b(j+2)-1, c_b(j)+s:c_b(j+1)-1+s) = H(1:r_b(j+2)-1, c_b(j):c_b(j+1)-1);
                
                K(1:r_b(j+2)-1, c_b(j):c_b(j+1)-1) = H(1:r_b(j+2)-1, c_b(j):c_b(j+1)-1);
                % Construct K.
                K(1:r_b(j+2)-1, c_b(j):c_b(j+1)-1) = nu(j-offset)*K(1:r_b(j+2)-1, c_b(j):c_b(j+1)-1) - ...
                    rho(j-offset)*[T(1:r_b(j+1)-1, c_b(j):c_b(j+1)-1); zeros(r,s)];
                % Construct H.
                H(1:r_b(j+2)-1, c_b(j):c_b(j+1)-1) = mu(j-offset)*H(1:r_b(j+2)-1, c_b(j):c_b(j+1)-1) - ...
                    eta(j-offset)*[T(1:r_b(j+1)-1, c_b(j):c_b(j+1)-1); zeros(r,s)];
                
                if param.balance ~= 0
                    % Balancing pencil
                    sd = diag(1./sqrt(sum(abs([ K(1:r_b(j+2)-1, c_b(j):c_b(j+1)-1); H(1:r_b(j+2)-1, c_b(j):c_b(j+1)-1) ]).^2)));
                    % Check case of zero column
                    sd(isinf(sd)) = 0;
                    K(1:r_b(j+2)-1, c_b(j):c_b(j+1)-1) = K(1:r_b(j+2)-1, c_b(j):c_b(j+1)-1)*sd;
                    H(1:r_b(j+2)-1, c_b(j):c_b(j+1)-1) = H(1:r_b(j+2)-1, c_b(j):c_b(j+1)-1)*sd;
                end
                
            else
                real_w = real(w);
                imag_w = imag(w);
                
                switch param.orth % for real part
                    case 'CGS'
                        hh = param.inner_product(real_w, V(:, 1:r_b(j+1)-1));
                        real_w = real_w - V(:, 1:r_b(j+1)-1)*hh;
                        H(1:r_b(j+1)-1, c_b(j):c_b(j)+s-1) = hh;
                        if param.reorth == 1
                            hh = param.inner_product(real_w, V(:, 1:r_b(j+1)-1));
                            real_w = real_w - V(:, 1:r_b(j+1)-1)*hh;
                            H(1:r_b(j+1)-1, c_b(j):c_b(j)+s-1) = H(1:r_b(j+1)-1, c_b(j):c_b(j)+s-1) + hh;
                        end
                        
                        [real_w, RR, E] = util_qr(real_w,param.inner_product); % Rank revealing qr.
                        Einv = E; Einv(E) = 1:length(E);
                        real_w = real_w(:,Einv);
                        RR = RR(Einv,Einv);
                        % 'CGS' end
                    case 'MGS'
                        for col = 1:j
                            hh = param.inner_product(real_w, V(:, r_b(col):r_b(col+1)-1));
                            real_w = real_w - V(:, r_b(col):r_b(col+1)-1)*hh;
                            H(r_b(col):r_b(col+1)-1, c_b(j):c_b(j)+s-1) = H(r_b(col):r_b(col+1)-1, c_b(j):c_b(j)+s-1) + hh;
                        end
                        if param.reorth == 1
                            for col = 1:j
                                hh = param.inner_product(real_w, V(:, r_b(col):r_b(col+1)-1));
                                real_w = real_w - V(:, r_b(col):r_b(col+1)-1)*hh;
                                H(r_b(col):r_b(col+1)-1, c_b(j):c_b(j)+s-1) = H(r_b(col):r_b(col+1)-1, c_b(j):c_b(j)+s-1) + hh;
                            end
                        end
                        
                        [real_w, RR, E] = util_qr(real_w, param.inner_product); % Rank revealing qr.
                        Einv = E; Einv(E) = 1:length(E);
                        real_w = real_w(:, Einv);
                        RR = RR(Einv, Einv);
                        % 'MGS' end
                end % switch param.orth
                
                % inf/nan check
                ssum = sum(sum([H(:,c_b(j):c_b(j)+s-1);RR]));
                if isinf(ssum) || isnan(ssum)
                    ME = MException('VerifyOutput:OutOfBounds', ...
                        ['Error: inf/nan occurred during orthogonalisation. Maybe pole ', ...
                       num2str(xi(j)), ' too close to eigenvalue of pencil?']);
                    throw(ME);
                end
                
                % Deflation test
                diagr = svd(RR);
                r = find(diagr > max(diagr(1)*param.deflation_tol, ...
                    norm(H(1:r_b(j+1)-1, c_b(j):c_b(j)+s-1),'fro')*param.deflation_tol), 1, 'last');
                if isempty(r), r = 0; end
                if r ~= s, bd1 = true; end
                
                r_b(j+2) = r_b(j+1) + r; c_b(j+1) = c_b(j) + s;
                RR(E(r+1:s),:) = [];
                H(r_b(j+1):r_b(j+2)-1,c_b(j):c_b(j)+s-1) = RR;
                real_w(:,E(r+1:s))=[];
                real_mem = E(r+1:s);
                V(:,r_b(j+1):r_b(j+2)-1) = real_w;
                inds = c_b(j)-1 + (1:s);
                param.column_deflation(inds) = ones(1,length(inds));
                j = j+1;
                
                switch param.orth % For imag part.
                    case 'CGS'
                        hh = param.inner_product(imag_w, V(:, 1:r_b(j+1)-1));
                        imag_w = imag_w - V(:, 1:r_b(j+1)-1)*hh;
                        H(1:r_b(j+1)-1, c_b(j):c_b(j)+s-1) = hh;
                        if param.reorth == 1
                            hh = param.inner_product(imag_w, V(:, 1:r_b(j+1)-1));
                            imag_w = imag_w - V(:, 1:r_b(j+1)-1)*hh;
                            H(1:r_b(j+1)-1, c_b(j):c_b(j)+s-1) = H(1:r_b(j+1)-1, c_b(j):c_b(j)+s-1) + hh;
                        end
                        
                        [imag_w, RR, E] = util_qr(imag_w,param.inner_product); % Rank revealing qr.
                        Einv = E; Einv(E) = 1:length(E);
                        imag_w = imag_w(:, Einv);
                        RR = RR(Einv, Einv);
                        % 'CGS' end
                        
                    case 'MGS'
                        for col = 1:j
                            hh = param.inner_product(imag_w, V(:, r_b(col):r_b(col+1)-1));
                            imag_w = imag_w - V(:, r_b(col):r_b(col+1)-1)*hh;
                            H(r_b(col):r_b(col+1)-1, c_b(j):c_b(j)+s-1) = H(r_b(col):r_b(col+1)-1, c_b(j):c_b(j)+s-1) + hh;
                        end
                        if param.reorth == 1
                            for col = 1:j
                                hh = param.inner_product(imag_w, V(:, r_b(col):r_b(col+1)-1));
                                imag_w = imag_w - V(:, r_b(col):r_b(col+1)-1)*hh;
                                H(r_b(col):r_b(col+1)-1, c_b(j):c_b(j)+s-1) = H(r_b(col):r_b(col+1)-1, c_b(j):c_b(j)+s-1) + hh;
                            end
                        end
                        
                        [imag_w, RR, E] = util_qr(imag_w,param.inner_product); % Rank revealing qr.
                        Einv = E; Einv(E) = 1:length(E);
                        imag_w = imag_w(:, Einv);
                        RR = RR(Einv, Einv);
                        % 'MGS' end
                end % switch param.orth
                
                % inf/nan check
                ssum = sum(sum([H(:,c_b(j):c_b(j)+s-1);RR]));
                if isinf(ssum) || isnan(ssum)
                    ME = MException('VerifyOutput:OutOfBounds', ...
                        ['Error: inf/nan occurred during orthogonalisation. Maybe pole ', ...
                       num2str(xi(j)), ' too close to eigenvalue of pencil?']);
                    throw(ME);
                end
                
                % Deflation test
                diagr = svd(RR);
                r = find(diagr > max(diagr(1)*param.deflation_tol, ...
                    norm(H(1:r_b(j+1)-1, c_b(j):c_b(j)+s-1),'fro')*param.deflation_tol), 1, 'last');
                if r ~= s, bd1 = true; end
                if isempty(r), r = 0; end
                
                r_b(j+2) = r_b(j+1) + r; c_b(j+1) = c_b(j) + s;
                RR(E(r+1:s),:) = [];
                H(r_b(j+1):r_b(j+2)-1,c_b(j):c_b(j)+s-1) = RR;
                imag_w(:,E(r+1:s)) = [];
                V(:,r_b(j+1):r_b(j+2)-1) = imag_w;
                inds = c_b(j)-1 + (1:s);
                param.column_deflation(inds) = ones(1,length(inds));
                
                % Deflation in real and complex part may not match.
                if not(isempty(E(r+1:s)))
                    if isempty(real_mem) || real_mem ~= E(r+1:s)
                        bd2 = true;
                        V(:,r_b(j):r_b(j+2)-1) = zeros(N, r_b(j+2)-r_b(j));
                        H(1:r_b(j+2)-1, c_b(j-1):c_b(j+1)-1) = zeros(r_b(j+2)-1,c_b(j+1)-c_b(j-1));
                        r_b(j+1) = r_b(j); r_b(j+2) = r_b(j);
                        c_b(j+1) = c_b(j-1); c_b(j) = c_b(j-1);
                        j = j+1;
                        continue;
                    end
                end
                
                % Setting up the decomposition.
                R(1:r_b(j+2)-1, c_b(j-1)+s:c_b(j+1)-1+s) = H(1:r_b(j+2)-1, c_b(j-1):c_b(j+1)-1);
                rmu = real(mu(j-1-offset)); imu = imag(mu(j-1-offset));
                d1 = c_b(j) - c_b(j-1);
                d2 = c_b(j+1) - c_b(j);
                
                % Construct 2x2 block permutation matrix.
                cmu = rmu*eye(2*s);
                cmu(1:s,s+1:2*s) = imu*eye(s);
                cmu(s+1:2*s,1:s) = -imu*eye(s);
                % Remove rows and columns linked to deflated.
                cmu(s+d2+1:2*s,:) = []; cmu(:,s+d2+1:2*s) = [];
                cmu(d1+1:s,:) = []; cmu(:,d1+1:s) = [];
                
                T1 = real(T(1:r_b(j)-1, c_b(j-1):c_b(j)-1));
                T2 = imag(T(1:r_b(j)-1, c_b(j-1):c_b(j)-1));
                
                % Construct modified continutation matrix.
                rcnt = [T1 T2; zeros(r_b(j+2)-r_b(j), d1+d2)];
                
                K(1:r_b(j+2)-1, c_b(j-1):c_b(j+1)-1) = H(1:r_b(j+2)-1, c_b(j-1):c_b(j+1)-1);
                K(1:r_b(j+2)-1, c_b(j-1):c_b(j+1)-1) = nu(j-offset)*K(1:r_b(j+2)-1, c_b(j-1):c_b(j+1)-1) - ...
                    rho(j-offset)*rcnt;
                H(1:r_b(j+2)-1, c_b(j-1):c_b(j+1)-1) = H(1:r_b(j+2)-1, c_b(j-1):c_b(j+1)-1)*cmu - ...
                    eta(j-offset)*rcnt;
                
                % Balancing pencil
                if param.balance ~= 0
                    sd = diag(1./sqrt(sum(abs([ K(1:r_b(j+2)-1, c_b(j-1):c_b(j+1)-1); H(1:r_b(j+2)-1, c_b(j-1):c_b(j+1)-1) ]).^2)));
                    % Check case of zero column
                    sd(isinf(sd)) = 0;
                    K(1:r_b(j+2)-1, c_b(j-1):c_b(j+1)-1) = K(1:r_b(j+2)-1, c_b(j-1):c_b(j+1)-1)*sd;
                    H(1:r_b(j+2)-1, c_b(j-1):c_b(j+1)-1) = H(1:r_b(j+2)-1, c_b(j-1):c_b(j+1)-1)*sd;
                end
                
            end % realopt
            
            if bd1 == true && bd2 == false
                columns_deflated = s-r;
                if columns_deflated == 1
                    warning('rat_krylov: 1 column deflated');
                elseif columns_deflated == s
                    warning('rat_krylov: all columns deflated');
                else
                    warning(['rat_krylov: ', num2str(s-r), ' columns deflated']);
                end
                bd1 = false;
            end
            if bd2 == true
                warning(['rat_krylov: deflation occurred during real variant, ' ...
                    'real and complex part do not match, poles ', num2str(xi(j-2)), ...
                    ' and ', num2str(xi(j-1)), ' removed.']);
                bd2 = false;
            end
            
            if param.waitbar, waitbar(j/m, hwb), end
            j = j+1;
        end % while j <= m
        
        if param.waitbar, close(hwb), end
        
        % When extendng, the previous block sizes are unknown so we replace
        % with NaN.
        if length(r_b) == length(xi)+2
            param.blocksizes = diff(r_b);
        else
            param.blocksizes = [NaN, diff(r_b(end - length(xi):end))];
        end
        
        % Remove excess zero rows and columns.
        V(:,r_b(j+1):size(V,2)) = [];
        K(r_b(j+1):size(K,1),:) = []; K(:,c_b(j):size(K,2)) = [];
        H(r_b(j+1):size(H,1),:) = []; H(:,c_b(j):size(H,2)) = [];
        %T(r_b(j+1):size(T,1),:) = []; T(:,c_b(j):size(T,2)) = [];
        R(r_b(j+1):size(R,1),:) = []; R(:,c_b(j)+s:size(R,2)) = [];
        param.column_deflation(size(H,2)+1:end) = [];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Nested functions of rat_krylov. %%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function continuation_pair
            j_loc = 1+floor((j-1)/param.p)*param.p;
            
            if ~isnan(param.continuation_root)
                if isinf(param.continuation_root)
                    rho(j-offset) = 0;
                    eta(j-offset) = 1;
                else
                    rho(j-offset) = 1;
                    eta(j-offset) = param.continuation_root;
                end
            end
            
            switch param.continuation
                case 'almost-last'
                    j_loc = max(1, j-param.p+1);
                    T(r_b(j_loc):r_b(j_loc+1)-1, c_b(j):c_b(j)+s-1) = eye(r_b(j_loc+1)-r_b(j_loc));
                    
                case 'first'
                    T(1:s, c_b(j):c_b(j)+s-1) = eye(s);
                    
                case 'last'
                    t = r_b(j_loc+1) - r_b(j_loc);
                    T(r_b(j_loc):r_b(j_loc+1)-1, c_b(j):c_b(j)+t-1) = eye(t);
                    
                case 'ruhe'
                    if j_loc == 1
                        T(r_b(j_loc):r_b(j_loc+1)-1, c_b(j):c_b(j)+s-1) = eye(s);
                    else
                        [Q, ~] = qr(K(1:r_b(j_loc+1)-1, param.column_deflation)*mu(j-offset) - ...
                            H(1:r_b(j_loc+1)-1, param.column_deflation)*nu(j-offset));
                        T(1:r_b(j_loc+1)-1, c_b(j):c_b(j)+s-1) = Q(:, end-s+1:end);
                    end
                    
                case 'near-optimal'
                    % Get auxiliary continuation vector. (Other choices may be considered.)
                    assert(s == 1 , 'near-optimal currently unavailable for block Krylov spaces')
                    if j_loc == 1
                        T(j_loc, j) = 1;
                    else
                        [Q, ~] = qr(K(1:j_loc, 1:j_loc-1)*mu(j-offset) - ...
                            H(1:j_loc, 1:j_loc-1)*nu(j-offset));
                        T(1:j_loc, j) = Q(:, end);
                    end
                    
                    % Approximate next vector.
                    vwhat = AB.multiply(rho(j-offset), eta(j-offset), V(:, 1:j_loc)*T(1:j_loc, c_b(j):c_b(j)+s-1));
                    
                    if param.continuation_bounds && param.p == 1
                        vshat = vwhat;
                    end
                    
                    vwhat = param.continuation_solve(AB, nu(j-offset), mu(j-offset), ...
                        vwhat, param);
                    
                    if param.continuation_bounds && param.p == 1
                        % Residual according to [3, eq. (3.8)].
                        vshat = AB.multiply(nu(j-offset), mu(j-offset), vwhat) - vshat;
                    end
                    
                    % CGS with reorthogonalization.
                    vchat  = param.inner_product(vwhat, V(:, 1:j_loc));
                    vvhat  = vwhat - V(:, 1:j_loc)*vchat;
                    
                    vchat2 = param.inner_product(vvhat, V(:, 1:j_loc));
                    vvhat  = vvhat - V(:, 1:j_loc)*vchat2;
                    vchat  = vchat + vchat2;
                    
                    vchat(end+1, 1) = param.induced_norm(vvhat);
                    vvhat = vvhat/vchat(end);
                    
                    % Form approximate pencil, cf. [3, eq. (3.4)].
                    K(1:j_loc+1, j) = nu(j-offset)*vchat - rho(j-offset)*T(1:j_loc+1, j);
                    H(1:j_loc+1, j) = mu(j-offset)*vchat - eta(j-offset)*T(1:j_loc+1, j);
                    
                    if param.continuation_bounds && param.p == 1
                        % Scaled error given by [3, eq. (3.12)].
                        vfhat = -AB.solve(nu(j-offset), mu(j-offset), vshat)/vchat(end);
                        % Corresponding norm, and norm of the local projection.
                        Fhat(j) = abs(param.induced_norm(vfhat));
                        fhat(j) = abs(param.induced_norm(...
                            [V(:, 1:j) vvhat] * ...
                            param.inner_product(vfhat, [V(:, 1:j) vvhat])));
                    end
                    
                    ind = 1;
                    [XX, DD] = eig(H(1:j_loc, [1:j_loc-1 j]), K(1:j_loc, [1:j_loc-1 j]));
                    
                    if isinf(DD(ind, ind))
                        rho(j-offset) = 0;
                        eta(j-offset) = 1;
                    else
                        rho(j-offset) = 1;
                        eta(j-offset) = DD(ind, ind);
                    end
                    
                    if isreal(mu(j-offset)) && ...
                            isreal(nu(j-offset)) && ...
                            isreal(H(1:j_loc, 1:j_loc-1)) && ...
                            isreal(K(1:j_loc, 1:j_loc-1)) && ...
                            ~isreal(eta(j-offset))
                        % Not looking for smallest imaginary part in relative
                        % sense to avoid zero.
                        [~, ind] = min(abs(imag(diag(DD))));
                        eta(j-offset) = DD(ind, ind);
                        eta(j-offset) = real(eta(j-offset));
                        XX = real(XX);
                    end
                    
                    % Scaling factor given by [3, eq. (3.6)]. (Can actually be excluded.)
                    gamma = XX(end, ind)*(rho(j-offset)*H(j_loc+1, j) - ...
                        eta(j-offset)*K(j_loc+1, j));
                    
                    LL = mu(j-offset)*K(1:j_loc, [1:j_loc-1 j]) - ...
                        nu(j-offset)*H(1:j_loc, [1:j_loc-1 j]);
                    
                    % Near-optimal continuation vector given by [3, eq. (3.7)].
                    T(1:j_loc, j) = LL*XX(:, ind);
                    T(1:j_loc, j) = T(1:j_loc, j)/gamma;
                    K(1:j+1, j) = 0;
                    H(1:j+1, j) = 0;
            end % switch param.continuation
        end % continuation_pair
    end % run_krylov

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function rerun_krylov
        r_b = param.r_b; c_b = param.c_b;% Retrieve block indexing.
        while j <= m
            % Computing the continuation combination.
            if isreal(xi(j)) || realopt == 0
                T(1:r_b(j+1)-1, c_b(j):c_b(j+1)-1) = K(1:r_b(j+1)-1, c_b(j):c_b(j+1)-1)*mu(j) - H(1:r_b(j+1)-1, c_b(j):c_b(j+1)-1)*nu(j);
                %U(1:j, j) = U(1:j, j)/(nu(j)*eta(j)-mu(j)*rho(j));
            elseif realopt == 1 || realopt == 2
                H1 = H(r_b(j+1):r_b(j+3)-1, c_b(j):c_b(j+2)-1);
                K1 = K(r_b(j+1):r_b(j+3)-1, c_b(j):c_b(j+2)-1);
                
                [HH, KK, QQ, ZZ] = qz(H1, K1); % HH=Q*H1*Z, KK=Q*K1*Z
                ee = ordeig(HH,KK);
                select = zeros(length(ee),1); select(imag(ee)<0) = 1;
                % Selectc orresponds to ordered poles in util_pencil_poles.
                [~,~,QS,ZS] = ordqz(HH,KK,QQ,ZZ,select); % KKK=QS*K1*ZS, HHH=QS*H1*ZS
                
                HHH = H(1:r_b(j+3)-1, c_b(j):c_b(j+2)-1)*ZS;
                HHH(r_b(j+1):r_b(j+3)-1, :) = QS*HHH(r_b(j+1):r_b(j+3)-1, :);
                KKK = K(1:r_b(j+3)-1, c_b(j):c_b(j+2)-1)*ZS;
                KKK(r_b(j+1):r_b(j+3)-1, :) = QS*KKK(r_b(j+1):r_b(j+3)-1, :);
                
                T(1:r_b(j+1)-1, c_b(j):c_b(j+1)-1) = KKK(1:r_b(j+1)-1,1:s)*mu(j) - HHH(1:r_b(j+1)-1,1:s)*nu(j);
                j = j+1;
                T(1:r_b(j+1)-1, c_b(j):c_b(j+1)-1) = KKK(1:r_b(j+1)-1,1+s:2*s)*mu(j) - HHH(1:r_b(j+1)-1,1+s:2*s)*nu(j);
                j = j-1;
            end
            
            % Next vector and orthogonalization.
            if isreal(xi(j)) || realopt == 0
                w = V(:, 1:r_b(j+1)-1)*T(1:r_b(j+1)-1, c_b(j):c_b(j+1)-1);
                w = AB.multiply(rho(j), eta(j), w);
                w = AB.solve(nu(j), mu(j), w);
                
                r = eta(j)*K(1:r_b(j+2)-1, c_b(j):c_b(j+1)-1) - rho(j)*H(1:r_b(j+2)-1, c_b(j):c_b(j+1)-1);
                %r = r/(nu(j)*eta(j)-mu(j)*rho(j))
                
                % MGS simulation.
                w = w - V(:, 1:r_b(j+1)-1)*r(1:r_b(j+1)-1,:);
                
                V(:, r_b(j+1):r_b(j+2)-1) = w / r(r_b(j+1):r_b(j+2)-1,:);
                
            else
                w = V(:, 1:r_b(j+1)-1)*T(1:r_b(j+1)-1, c_b(j):c_b(j+1)-1);
                w = AB.multiply(rho(j), eta(j), w);
                w = AB.solve(nu(j), mu(j), w);
                
                r = eta(j)*KKK(1:r_b(j+2)-1, 1:s) - rho(j)*HHH(1:r_b(j+2)-1,1:s);
                
                % MGS simulation.
                w = w - V(:, 1:r_b(j+1)-1)*r(1:r_b(j+1)-1, :);
                
                V(:, r_b(j+1):r_b(j+2)-1)= w / r(r_b(j+1):r_b(j+2)-1,:);
                
                j = j+1;
                
                w = V(:, 1:r_b(j+1)-1)*T(1:r_b(j+1)-1, c_b(j):c_b(j+1)-1);
                w = AB.multiply(rho(j), eta(j), w);
                w = AB.solve(nu(j), mu(j), w);
                
                r = eta(j)*KKK(1:r_b(j+2)-1, s+1:2*s) - rho(j)*HHH(1:r_b(j+2)-1, s+1:2*s);
                
                % MGS simulation.
                w = w - V(:, 1:r_b(j+1)-1)*r(1:r_b(j+1)-1, :);
                
                V(:, r_b(j+1):r_b(j+2)-1)= w / r(r_b(j+1):r_b(j+2)-1,:);
                
                V(:, r_b(j):r_b(j+2)-1) = V(:, r_b(j):r_b(j+2)-1)*QS;
                
                if realopt == 1, V(:, r_b(j):r_b(j+2)-1) = real(V(:, r_b(j):r_b(j+2)-1)); end
            end
            j = j+1;
        end % while j <= m
    end % rerun_krylov

warning on backtrace
end % rat_krylov

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [AB, V, K, H, N, j, m, ...
    s, rerun, param] = parse_argin(varargin)
%PARSE_ARGIN    Process the input argument list to rat_krylov.

msg = 'rat_krylov:parse_argin: ';
assert(nargin >= 3, 'more parameters needed');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Defines AB, N, and j. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check for parameters given
if isstruct(varargin{end})
    param_given = 1; % Check structure at end
elseif ischar(varargin{end})
    param_given = 1; % Check parameter at end
else
    param_given = 0;
end

inp = nargin-param_given;

if isstruct(varargin{1}) % Check pencil structure AB
    AB = varargin{1}; % AB index
    V_id = 2; % must be b or V
else
    if inp == 3
        V_id = 2; B_id = 0;
    elseif inp == 4
        if size(varargin{3}) == size(varargin{4})
            V_id = 2; B_id = 0;
        else
            V_id = 3; B_id = 2;
        end
        
    elseif inp == 5
        if size(varargin{4}) == size(varargin{5})
            V_id = 3; B_id = 2;
        else
            V_id = 2; B_id = 0;
        end
    elseif inp == 6
        V_id = 3; B_id = 2;
    else
        warning('parameter inputs invalid')
        return
    end
end

[N, j] = size(varargin{V_id}); % dimensions of our b or V

% No Structure
if ~isstruct(varargin{1}) % Check A is square
    assert(size(varargin{1}, 1) == N, [msg 'A has to be N-by-N where b has N rows'])
    assert(size(varargin{1}, 2) == N, [msg 'A has to be N-by-N where b has N rows'])
    if B_id == 0 % Create function handles for just matrix
        AB.multiply = @(rho, eta, x) rho*(varargin{1}*x) - eta*x;
        AB.solve    = @(nu, mu, x) (nu*varargin{1}-mu*speye(N))\x;
        AB.isreal = isreal(varargin{1});
    else % Check B is square, define function handles
        assert(size(varargin{2}, 1) == N, [msg 'B has to be N-by-N'])
        assert(size(varargin{2}, 2) == N, [msg 'B has to be N-by-N'])
        AB.multiply = @(rho, eta, x) ...
            rho*(varargin{1}*x) - eta*(varargin{2}*x);
        AB.solve    = @(nu, mu, x) ...
            (nu*varargin{1}-mu*varargin{2})\x;
        AB.isreal = isreal(varargin{1}) && isreal(varargin{2});
    end
end

if ~isfield(AB, 'isreal'), AB.isreal = false; end % Check parameter input for real version of Arnoldi.

% Structure
if ~all(isfield(AB, {'multiply', 'solve'})) % Check multiply and solve have been definied.
    assert(any(isfield(AB, {'A', 'B'})), ...
        [msg 'the pencil is null'])
    if ~isfield(AB, {'B'}) % If given just A
        assert(size(AB.A, 1) == N, [msg 'A has to be N-by-N'])
        assert(size(AB.A, 2) == N, [msg 'A has to be N-by-N'])
        AB.isreal = isreal(AB.A);
        AB.multiply = @(rho, eta, x) rho*(AB.A*x) - eta*x;
        AB.solve    = @(nu, mu, x) (nu*AB.A-mu*speye(N))\x;
    elseif ~isfield(AB, {'A'}) % If given just B
        assert(size(AB.B, 1) == N, [msg 'B has to be N-by-N'])
        assert(size(AB.B, 2) == N, [msg 'B has to be N-by-N'])
        AB.isreal = isreal(AB.B);
        AB.multiply = @(rho, eta, x) rho*x - eta*(AB.B*x);
        AB.solve    = @(nu, mu, x) (nu*speye(N)-mu*AB.B)\x;
    else % Given both
        assert(size(AB.A, 1) == N, [msg 'A has to be N-by-N'])
        assert(size(AB.A, 2) == N, [msg 'A has to be N-by-N'])
        assert(size(AB.B, 1) == N, [msg 'B has to be N-by-N'])
        assert(size(AB.B, 2) == N, [msg 'B has to be N-by-N'])
        AB.isreal = isreal(AB.A) && isreal(AB.B);
        AB.multiply = @(rho, eta, x) rho*(AB.A*x) - eta*(AB.B*x);
        AB.solve    = @(nu, mu, x)(nu*AB.A-mu*AB.B)\x;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Defines V, K, H, xi, m, rerun and param. %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rat_krylov_mode = nargin - V_id-param_given; % Check with functionality I,II or III.
rerun = 0;
r_b = zeros(1,0); c_b = zeros(1,0);


if ~param_given || ischar(varargin{end})
    param = struct;
    if strcmp('real', varargin{end}), param.real = 1; end
else
    param = varargin{end};
end

if isfield(param, 'orth')
    if ~strcmp(param.orth, 'CGS') && ~strcmp(param.orth, 'MGS')
        warning('Warning: param.orth invalid. Setting to default.')
        param.orth = 'MGS';
    end
else
    param.orth = 'MGS';
end

switch rat_krylov_mode
    case {1} % Rational Krylov.
        s = j;
        assert(rank(varargin{V_id})== size(varargin{V_id},2), 'b must have full column rank')
        j = 1;
        xi = varargin{V_id+1};
        if size(xi,1) == 0
            xi = zeros(1,0);
        else
            assert(size(xi, 1) == 1, ...
            [msg 'poles must be row vector'])
        end
        m  = length(xi);
        V = zeros(N, s*(m+1))*varargin{V_id}(1,1); % Support mp
        V(:, 1:s) = varargin{V_id};
        K = zeros(s*(m+1), s*(m))*varargin{V_id}(1,1); % Support mp
        H = zeros(s*(m+1), s*(m))*varargin{V_id}(1,1); % Support mp
        param.column_deflation = false(1,s*m);
        
    case {2} % Rerunning.
        rerun = 1;
        K = varargin{V_id+1};
        H = varargin{V_id+2};
        s = size(varargin{V_id},2);
        assert(size(K, 1) == size(H, 1), ...
            [msg 'K and H need to be of equal size'])
        assert(size(K, 2) == size(H, 2), ...
            [msg 'K and H need to be of equal size'])
        %assert((0 < sm) && (sm <= s), ...
        %[msg 'K has to be of size (m+sm)-by-m where 0 < sm <= s0'])
        V = zeros(N, size(K, 1));
        
        % Allow for multiple precision rerunning.
        if isa(K, 'sym') || isa(H, 'sym') || isa(varargin{V_id}, 'sym')
            V = sym(V);
        elseif isa(K, 'mp') || isa(H, 'mp') || isa(varargin{V_id}, 'mp')
            V = mp(V);
        end
        
        V(:,1:s) = varargin{V_id};
        [xi, r_b] = util_pencil_poles(K, H, s);
        if size(xi, 1) == 0
            xi = zeros(1,0);
        end
        c_b = 1:s:size(H,2)+1;
        m = length(xi);
        j = 1;
        
    case {3} % Extending.
        if ~isfield(param, 'extend')
            param.extend = 1;
        end
        s = param.extend; % Number of vectors we want to extend by.
        assert( s <= j, ...
            [msg 'param.extend is too large.'])
        xi = varargin{V_id+3};
        if size(xi, 1) == 0
            xi = zeros(1,0);
        end
        m  = length(xi);
        assert(size(varargin{V_id+1}, 1) == j, ...
            [msg 'K must have j rows'])
        assert(size(varargin{V_id+2}, 1) == j, ...
            [msg 'H must have j rows'])
        assert(size(varargin{V_id+2}, 2) == size(varargin{V_id+1}, 2), ...
            [msg 'K and H need to be of equal size'])
        
        V = zeros(N, (m*s)+j);
        [Hr, Hc] = size(varargin{V_id+1});
        K = zeros((m*s)+Hr, (m*s)+Hc);
        H = zeros((m*s)+Hr, (m*s)+Hc);
        V(:, 1:j) = varargin{V_id};
        K(1:Hr, 1:Hc) = varargin{V_id+1};
        H(1:Hr, 1:Hc) = varargin{V_id+2};
        
        param.column_deflation = logical([ones(1,Hc), zeros(1,m*s)]);
        
        if strcmp(param.orth, 'CGS')
            r_b = [1, Hr+1-s, Hr+1];
            c_b = [1, Hc+1];
        else
            if Hr > Hc
                r_b = [1:Hc, Hr+1-s,Hr+1];
                c_b = 1:Hc+1;
            else
                r_b = [1:Hr+1-s, Hr+1];
                c_b = [1:Hr-s, Hc+1];
            end
        end
        
        j = length(r_b)-1;
        m = m + j - 1;
        
    otherwise
        error([msg 'invalid number of input arguments']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isfield(param, 'reorth')
    if param.reorth ~= 0 && param.reorth ~= 1
        warning('Warning: param.reorth invalid. Setting to default.')
        param.reorth = 1;
    end
else
    param.reorth = 1;
end

if isfield(param, 'real')
    if param.real ~= 1 && param.real ~= 0
        warning('Warning: param.real invalid. Setting to default.')
        param.real = 0;
    end
else
    param.real = 0;
end

if isfield(param, 'refinement')
    if param.refinement ~= 0 && param.refinement ~= 1
        warning('Warning: param.refinement invalid. Setting to default.')
        param.refinement = 0;
    end
else
    param.refinement = 0;
end

if rat_krylov_mode == 2 && size(H, 2) > 1 && any(diag(H, -2*s))
    param.real = 1;
end

if param.real && (~AB.isreal || ~isreal(V(:, 1:s)))
    if rat_krylov_mode == 2
        param.real = 2;
    else
        param.real = 0;
        warning(['Warning: Ignoring ''real'' option with complex' ...
            ' matrices/vectors.']);
    end
end
if param.real && ~canonical_cplx(xi)
    param.real = 0;
    warning(['Warning: Ignoring ''real'' option as conj(poles(j))' ...
        ' == poles(j+1) failed']);
end

param.xi = xi;
param.moebius = poles_to_moebius(xi);
param.r_b = r_b; param.c_b = c_b;

if ~isfield(param, 'inner_product')
    param.inner_product = @(x, y) y'*x;
end
param.induced_norm  = @(x) sqrt(param.inner_product(x, x));

if isfield(param, 'waitbar')
    if param.waitbar ~= 1 && param.waitbar ~= 0
        warning('Warning: param.waitbar invalid. Setting to default.')
        param.waitbar = 0;
    end
else
    param.waitbar = 0;
end

if ~isfield(param, 'continuation')
    param.continuation = 'ruhe';
end

if ~isfield(param, 'continuation_root')
    param.continuation_root = NaN;
end

if ~isfield(param, 'continuation_solve')
    param.continuation_solve = @(AB, nu, mu, x, param) AB.solve(nu, mu, x);
end

if ~isfield(param, 'continuation_bounds')
    param.continuation_bounds = 0;
end

if ~isfield(param, 'continuation_m')
    param.continuation_m = 4;
end

if ~isfield(param, 'deflation_tol')
    param.deflation_tol = eps(1);
end

if isfield(param, 'balance')
    if param.balance ~= 0 && param.balance ~=1
        warning('Warning: param.balance invalid. Setting to default.')
        param.balance = 1;
    end
else
    param.balance = 1;
end

% To simulate parallel execution.
if ~isfield(param, 'p') % Set to default 1.
    param.p = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

assert(rerun || m < N, [msg 'number of poles cannot be greater' ...
    ' or equal to the size of the problem']);
assert(size(xi, 1) == 1 || size(xi,1) ==0 , ...
    [msg 'xi must be a 1-by-m row vector, or 0x0 double ']);
% assert(all(~any(tril(K(3:end, 1:end-1)))), ...
%         [msg 'H is not upper-Hessenberg']);
% assert(all(~any(tril(H(4:end, 1:end-1)))), ...
%         [msg 'K is not quasi-upper-Hessenberg']);
end % parse_argin


function y = canonical_cplx(xi)
% CANONICAL_CPLX Check if the poles xi are ordered canonically.

m = length(xi);
y = 1;
j = 1;
while j <= m
    if isreal(xi(j)) || isinf(double(xi(j)))
        j = j+1;
    else
        if j == m || (j < m && xi(j+1) ~= conj(xi(j)))
            y = 0;
        end % if
        j = j+2;
    end % if
end % while
end


function moebius = poles_to_moebius(xi)
% POLES_TO_MOEBIUS Moebius transformation with poles xi.
%
% Finite xi is replaced with (nu, mu) := (1, xi) and
% (rho, eta) := (0, 1),  and xi = inf is replaced by
% (nu, mu) := (0, 1) and (rho, eta) := (1, 0).

nu  = ones(1, length(xi));
mu  = xi;
rho = zeros(1, length(xi));
eta = ones(1, length(xi));

%nu(isinf(xi))  = 0; % does not work with vpa
nu(abs(xi)==inf) = 0;
%mu(isinf(xi))  = 1; % does not work with vpa
mu(abs(xi)==inf) = 1;
rho(abs(xi)>1) = 1;
eta(abs(xi)>1) = 0;

moebius = [nu; mu; rho; eta];
end