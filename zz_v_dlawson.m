function  [Rfun, relgap, dk, err, pol, zers, bestW] ...
    = zz_v_dlawson(x1, bigF, nd, deg_vec, maxit )
%
%  --------------------- vec_d_lawson --------------------
%         produce rational approximation of vector-valued function
%
%         min  MAX_j     ( |f1-p1/q|^2 + ... + |fk - pk/q|.^2 )_j
%              1<j<m
%
%         by solving
%               d2w = min w*( ||f1*q - p1||^2 + ... + ||fk*q - pk||^2 )
%               update w by lawson's method
% ---------------------------------------------------------
%
% Input:
%
%           x1                 sample points with length m
%           bigF           (m-by-n2) matrix [f_1,f_2,...,f_k]
%           nd               the degree of  denominator polynomial
%           deg_vec        degrees (n1,n2,...,nk) of numerators  p1,...,pk
%           maxit          maximun lawson's iteration step
% Output:
%
%          Rfun             cellarray with Rfun{1} = r1, ... , Rfun{ell} = r_ell
%                               each r_k = p_k/q ~ f_k is
%                           function handle with  degree (nk,nd)
%          pol              zeros of q
%          zers             cellarray with zers{k} is zeros of pk
%          dk               the sequence of the dual function values d2w
%          err              the record of maximun F-norm errors, i.e.,
%                               err(k) = MAX_j ||bigF(j,:)-bigR(j,:)||_{fro};
%          relgap           relative daulity gap at the computed approximant;
%                           generally, relgap < 0.01 means the computed can be viewed
%                           as the minimax approximant;
%
%
% remark: this code is also suitable for the scalar-value function, as well
%         as the degree-specified duplexer problem
% ------------------------------------------------------------

% Reference:
% [1] L.-H. Zhang, Y.-N. Zhang, C. Zhang and S. Han, Rational minimax 
%     approximation of matrix-valued functions, 
%     2025, URL https://arxiv.org/pdf/2508.06378v2.
% [2] L.-H. Zhang and C. Zhang, Rational minimax approximations for 
%     matrix-valued functions:Existence, optimality and algorithms, 
%     URL https://arxiv.org/pdf/2607.22576.
% ----------------------------------------------------------------

x1 = x1(:);     m = length(x1);

n2 = size(bigF,2);                   % bigF must be m-by-n2 matrix

nn = length(deg_vec);
if nn ~= n2
    error('Input #of F and #of degree do not matching!')
end


maxit = maxit + 1; 
% -------------  N largest degree of polynomial by Arnoldi

N = max([deg_vec(:);  nd]);

% initial weight
w = ones(m,1);                v = w./norm(w,2);

err = [];      dk = [];

tol = 0.01;


 
% -------------------- lawson's iteration ---------------------------


% maxit = maxit + 1;       
err0 = inf;

for kk = 1:maxit
    
 
    %   step 1.  ------------------------------------------
    %            basis functions  Q // span(Q) = krylov(diag(X),v);
    
    [Q,H] = zz_arnoldi(x1, v, N+1) ;
    
    Qd = Q(:,1:nd+1);           % basis matrix of q
    
    
    %   step 2. ------------------------------------------
    %           compute d(w)   ( vector-value function )
    %
    
    %         PFQ = [ (I-Q1*Q1')*(f1.*Qd);
    %                 (I-Q2*Q2')*(f2.*Qd); ...
    %                 (I-Qk*Qk')*(fk.*Qd)  ]
    
    PFQ = zeros(m*n2, nd+1);
    for j = 1:n2
        
        nj = deg_vec(j);
        Qj = Q(:,1:nj+1);
        
        fQ = bigF(:,j).*Qd;
        
        mj = fQ - Qj*(Qj'*fQ );
        
        II = (j-1)*m + (1:m);
        PFQ(II,1:nd+1) =  mj;
        
    end
    
    
%     [~,ss,cc] = svds(PFQ,1,'smallest');   %    TLS by svd
%    -------  more accurate but slowly, 
                [~,ss,cc] = svd(PFQ,0);     ss = ss(end,end);
%     cond(PFQ)
    
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    % % %                    niubi porjection (on the way ? )
    
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    
    
    %   step 3. ------------------------------------------
    %                 compute q and p1,p2,... with coeffs
    %                         b,and Ak{1,2,...}
    
    b = cc(:,end);
    %     %   ----------- constraint ----------- scaled b s.t. |w*q| = 1
 
%     cst =  norm(abs(Qd*b)),
%     b = b./cst;
    qq = Qd*b;
    
    %   --------------------------- q = Q*b  // pj = Q*a{j}
    Ak = cell(1,n2);       bigR = zeros(m,n2);          bigP = zeros(m,n2);
    for j = 1:n2
        nj = deg_vec(j);
        Qj = Q(:,1:nj+1);
%         aj = Qj'*(bigF(:,j).*qq);
        aj = Qj\(bigF(:,j).*qq);
        
        Ak{j} = aj(:);
        
        pj = Qj*aj;
        bigR(:,j) = pj./qq ;
        bigP(:,j) = pj ;
    end
    
%     
%     II = find(abs(qq)<1e-14) , 
%     
%     bigR(II,:) = bigF(II,:); 
%     
    
    
     
    
    d2w = ss ;      dk = [dk, d2w];
    
    %  step 4.  ---------------- compute error ------------------
    %                 rr =  MAX_j \| bigF(j,:) - bigR(j,:) \|_F,
    
    res = bigF - bigR ;
    rr =  sqrt( sum( abs(res).^2, 2) ) ;
    eta = norm(rr(:), 'inf');
    err = [err,  eta ];
    
    if eta < err0
        bestB = b;   bestA = Ak;
        bestW = w;   bestK = kk;
    end
    
    
    
    %   step 5. ------- check error and relative gap -------------
    
    relgap = abs(err(end)- dk(end))/err(end);
    
    if relgap < tol
        fprintf('\n vec_d_lawson stop at step %d with relative gap is %.4f \n ', kk, relgap);
        bestB = b;   bestA = Ak;
        bestW = w;   bestK = kk;
        break,
    end
    
    
    if  (err(end) < 2e-15)
        bestB = b;   bestA = Ak;
        bestW = w;   bestK = kk;
        fprintf('\n max error is %.4e at step k = %d \n', err(end), kk );
        break,
    end



    
    
    %   step 6. ----------- update weight ----------------------
    %                    lawson's method, w = w_old.*rr
    
    
    w = w.*abs( rr );  w = abs(w)./norm(w,1);    v = w.^0.5;
    
    %         w = w./abs(qq);   w = w./norm(w,1);   v = w.^0.5 ;   % ssk
    
    
    
end

 



%    ---------------------output ---------------------------
%           function handle, poles and zeros

Ak = bestA;   b = bestB;    v = bestW.^0.5;

[~,H] = zz_arnoldi(x1, v, N+1) ;

dk = dk(1:bestK);  err = err(1:bestK);


Rfun = cell(1, n2);
for j = 1:n2
    aj = Ak{j};
    Rfun{j} = @(zz) zz_ratevalA(zz, aj, b, H) ;
end


if nd > 0

[pol, ~] = zz_roots(H, b);

else
    pol = [];
end

zers = cell(1, n2);
for j = 1:n2
    aj = Ak{j}; 
    [xx,~] = zz_roots(H, aj);
    zers{j}  = xx;
end



end







% %  --------------------- subroutine ---------------------


function [Q,H,k,v_norm] = zz_arnoldi(A,v,n)
% Arnoldi process with the modified Gram-Schmidt to generate A,Q,H.
% Re-orthogonalization was used for illcondtioned matrix.
% Full Orthogonalization Method: A*Q = Q*H
%  
%-----------------------------------------------------------
%  k-step Arnoldi:    A Q(:,1:k) = Q H, Q(:,1) = v0/norm(v0) 
%-----------------------------------------------------------
%
%  Input:
%
%      A      (n-by-n)  matrix  or  n-vector
%      v       n-vector,  
%      n      # of Arnoldi steps requested
%  Output:
%      Q       n-by (k+1)  matrix  that Q'*Q = I 
%      H       Hessenberg matrix  (k+1)\times k
%      k      # of Arnoldi steps   act

opts.mGS = 1;            % modified GS process
opts.reorth = 1;       % 1 = double orthogonal
opts.select = 1;       % selective reorth when reorth = 0 and select = 1;

if isvector(A)
    m = length(A); 
    A = sparse(1:m,1:m, A(:) ); 
end



m = size(A,1);    

if isempty(v)
    v = ones(m,1);
end



Q = zeros(m,n);    H = zeros(n+1,n);
v_norm = norm(v,2);   v = v./v_norm;
Q(:,1) = v;

for k = 1:n
    w = A*Q(:,k);
    if opts.mGS == 0
        %    basic gram schmidt
        H(1:k,k) = Q(:,1:k)'*w;    % projection coefs
        w = w - Q(:,1:k)*H(1:k,k); % w orth to Q(:,1:k)
    else
        %            modified Gram-Schmidt
        for i = 1:k
            v = Q(:,i);
            H(i,k) = v'*w;
            w = w - v*H(i,k);
        end
    end
    %   re-orth
    if  opts.reorth == 1
        for i = 1:k
            tmp = Q(:,i)'*w;
            w = w - Q(:,i) * tmp;
            H(i,k) = H(i,k) + tmp;
        end

    elseif opts.select == 1

        for i = 1:k
            tmp = Q(:,i)'*w;
            if abs(tmp)> 2*eps
                w = w - Q(:,i) * tmp;
                H(i,k) = H(i,k) + tmp;
            end
        end
    end

    H(k+1,k) = norm(w,2);
    Q(:,k+1) = w / H(k+1,k);

    if H(k+1,k) < 2e-15
        fprintf('\n arnoldi process break down at step %d \n',k); 
        H = H(1:k,1:k); 
        Q = Q(:,1:k); 
        break
    end

end
 




end








