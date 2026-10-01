function [Rfun, mJ,  J,  pol,  zers,  err ] = ...
    zz_v_aaa_lawson(t, bigF, nd, deg_vec, maxit)
%
%   --------------- vec_AAA_lawson -------------------------
%    rational approximation : (no objective)
%               min | f1-p1/q | + |f2-p2/q| + ...+ |fk - pk/q| 
% 
%   ---------------------------------------------------------
% Input:
%
%           t                 sample points with length M
%           bigF           (M by ell) matrix [f_1,f_2,...,f_ell], each f is M-vector
%           nd               the degree of  denominator polynomial
%           deg_vec      degrees (n1,n2,...,n_ell) of numerators  p1,...,p
%           maxit          maximun lawson's iteration step
% Output:
%
%          Rfun             cellarray with Rfun{1} = r1, ... , Rfun{ell} = r_ell
%                              each r_k = p_k/q ~ f_k is function handle with  degree (n,n)
%          pol               zeros of q
%          zers             cellarray with zers{k} is zeros of p_k
%          mJ                 index of support point in sample pts t
%          J                   (1:M)\mJ
%
% ---------------------------------------------------------------------
% remark:
% (1) the method is different from zz_v_aaa when  nd~= deg_vec(n1,n2,...)
% (2) if nd = n1 = n2 = ..., it returns the results of
%                                          v_aaa(t,bigF,nd,'lawson',maxit)


% Applicaitons: rational approx for s11 = N/D  S21 = Pt/ D; S31 = Pr / D;
%                         where D is polynomial(degree nd),
%                         N,Pt,Pr are polynomials with degrees deg_vec
%                         (n1,n2,n3) given and fixed
% ---------------------------------------------------------------------
% References:
% [1] Zhang Lei-Hong, Ya-Nan Zhang, Linyi Yang, Ruwu Xiao, A AAA-type
%     algorithm for the microwave duplexer filtering, 2024.
% [2] J.-P. Berrut and H. D. Mittelmann.
%       Matrices for the direct determination of the barycentric weights
%       of rational interpolation. J. Comput. Appl. Math., 78(2):355??370, 1997.
% ---------------------------------------------------------------------

Z = t(:);      M = length(Z);

[MM,~] = size(bigF);
if MM~=M
    error('\n size of the input data (t, bigF) not match! \n ')
end


% step 0.  ---------  get support pts by v_AAA  ----------------

if  max(deg_vec) > nd
    error('\n degree of numerator must be less than that of denominator ! \n ')
elseif min(deg_vec) == nd
    [ Rfun,  mJ,  J,   pol,  zers,  err] = zz_v_aaa(t, bigF, nd, 'lawson',maxit );
    return
else
    [ Rfun0,  mJ,  J,   pol0,  zers0,  err0] = zz_v_aaa(t, bigF, nd );
end



%    check if more support pts (larger degree nd) input,  this code works
%    in the new degree nd, and deg_vec

if err0 < 1e-14
    nd = length(mJ) - 1;    deg_vec = min(nd,  deg_vec );
    fprintf('\n  the actually compute degrees nd =  %d  and deg_vec = \n ', nd )
    disp( deg_vec)
    if min(deg_vec) == nd
        [ Rfun,  mJ,  J,   pol,  zers,  err] = zz_v_aaa(t, bigF, nd, 'lawson',maxit );
        return
    end
end


%    step 1.-------------- delete support pts form sample pts -------------- 
%    
% % Z is support pts
Z = t(mJ);

% lZ = length(Z);

 
t = t(J);
 
bigF =   bigF(J,:) ;

[M, ell] = size(bigF);

% cauchy matrix
C = 1./ (t(:) - Z.');

bigM = [];
bigr = [];
QQ = cell(1,ell);
for ll = 1:ell
    nl = deg_vec(ll);
    mj = getbasis(Z,nd,nl);
    QQ{ll} = mj;
    bigM = blkdiag(bigM,C*mj);
    bigr = vertcat(bigr, bigF(:,ll).*C );
end
bigM0 = [bigM, -bigr];

%   lawson iteration initial weight 

bigS = 1;       wt = cell(1,ell);
for ll = 1:ell
    wt{ll} = 1;
end

bigR = zeros(M,ell);


% alpA = [];   betA = []; 
% % step 2. ----------------- lawson's iteration  -----------------
maxit = maxit + 1;       err0 = inf; 

for kk = 1 : maxit
    %   bigS:  iteration matrix (weight matrix W)
 
            bigM = bigS*bigM0;
            % Total least square by "economy SVD "
%             [~,~,V] = svd(bigM,0);    v = V(:,end);
                [~,~,v] = svds(bigM,1,'smallest');
            %     v: minmium sigular vector

            %  coeffs of numerator
            alp = cell(1,ell);
            for ll = 1:ell
                Ql = QQ{ll};
                nl = deg_vec(ll);
                mj = v(1:nl+1);
                alp{ll} = Ql*mj;
                v(1:nl+1) = [];
            end
            %  coeffs of denominator
            bet = v;

            % approximation of bigR ~= bigF
            bigD = C*bet;
            for ll = 1:ell
                mj = alp{ll};
                bigR(:,ll) = (C*mj)./bigD;
            end


            e1  = abs(bigF - bigR);

            % check error and tol
            errA = norm(e1(:),"inf");
            if errA < 4e-15
                alpA = alp;    betA = bet ;
                fprintf('\n lawson iteraton achives max error = %.4e at step %d \n', errA, kk)
                break
            end

            %     remain the best "alp and bet"
            if errA < err0
                alpA = alp;    betA = bet ;
                err0  = errA;   
            end

            %             % compute error and update weight
            %             for ll = 1:ell
            %                 err(:,ll) = err(:,ll)./ norm( err(:,ll), "inf") ;
            %             end

            www = zeros(M,ell);
            for ll = 1:ell
                mj = wt{ll};
                mj = mj.*e1(:,ll) ;
                %                 mj = mj./norm(mj,1);
                wt{ll} = mj;
                www(:,ll) = mj;
            end

            cst = norm(www(:),1);

            www = www(:)./cst;
            for l = 1:ell
                wt{l} = wt{l}./cst;
            end

            % sparse weight matrix
            bigS = sparse(1:ell*M,1:ell*M, sqrt( www(:) )  );

 
end

 
% step 3 -----------------------------------------------------------
 %         output rational function and poles and zeros for lawson's result

    err  =  errA; 
    %  function handle
    Rfun = cell(1,ell);

    wj = betA;    zj = Z;

    for ll = 1:ell
        fj1 = alpA{ll}  ;  
        fj1 = fj1 ./ wj ; 
        R1 = @(zz) zz_ratevalB(zz, zj, fj1, wj);
        Rfun{ll} = R1;
    end


    %  poles and zeros


    zers = cell(1,ell);
    for ll = 1:ell
        fj1 = alpA{ll}; 
         fj1 = fj1 ./ wj ; 
        [pol,~, mj ] = zz_rootsB(zj,fj1, wj );
        zers{ll} = mj;
    end

 
end











% %  --------------------- subroutine ---------------------
% %   get Q and make C*Q*alp = sum_j=0^n (Q*alp)_j / (z - z_j)
% %   form the zm degree polynomial
%
function Q = getbasis(Z,n,m)
% Z : support pts
% R(n,n) \to R(m,n)
if n==m
    Q = eye(n+1);
else
    V = Z.^(0:n-m-1);
    [Q,~] = qr(V);
    Q = Q(:,n-m+1:end) ;
end

end






