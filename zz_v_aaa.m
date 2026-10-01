function [ Rfun,  mJ,  J,   pol,  zers,  err] = zz_v_aaa(t, bigF, n, varargin )
%
% zz_v_aaa is a AAA-type method with Lawson's iteration;
% it produces vector-valued  rational approximants (n,n)type
% as well as the corresponding poles and zeros
%      --------------- vec_AAA -------------------------
%     rational approximation : 
%               min | f1-p1/q | + |f2-p2/q| + ...+ |fk - pk/q| 
%
% the modified code 'zz_v_aaa_lawson' produce rational approximation
% (m1,n), ..., (m_ell,n) type for user-specific.
%
% Input:
%
%           t                 sample points with length M
%           bigF           (M by ell) matrix [f_1,f_2,...,f_ell], each f is M-vector
%           n                 the degree of numerator/denominator polynomial
%
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
% References:
% [1] Zhang Lei-Hong, Ya-Nan Zhang, Linyi Yang, Ruwu Xiao, A AAA-type
%     algorithm for the microwave duplexer filtering, 2024.
% [2] Y. Nakatsukasa, O. Sete, and L. N. Trefethen. The AAA algorithm for
%     rational approximation. SIAM J. Sci. Comput., 40:A1494-A1522, 2018.
% [3] P. Lietaert, K. Meerbergen, J. Perez, and B. Vandereycken. Automatic
%     rational approximation and linearization of nonlinear eigenvalue
%     problems. IMA J. Numer. Anal., 42:1087-1115, 2022.
% ---------------------------------------------------------------------


Z = t(:);      M = length(Z);

[MM,ell] = size(bigF);
if MM~=M
    error('size of the input data (t, bigF) not match! ')
end

J = 1:M;

% inital support pts
z = [];
% smallF = [f1,f2,f3,... ]

smallF = [];

% f1j = [];  f2j = []; f3j = [];

% inital cauchy matrix
C = [];
bigR = mean(bigF);
bigR = repmat(bigR, M,1);

for k = 1 : n + 1
    %     for scalar fun
    %     [~,j] = max(abs( F - R ) );
    %     choose support pts by greedy algorithm

    % L2 norm (F norm)
    [~,j] = max(sum(abs(bigF - bigR).^2, 2 ) );
    %           L1 norm
    %               [~,j] = max(  sum(abs(bigF - bigR ), 2)   );
    %          Inf norm
    %               [mj,~] = max( max(abs(bigF - bigR )) ) ;
    %               [j,~] = find(abs(bigF - bigR) == mj);

    z = [z; Z(j)];
    smallF =  vertcat( smallF, bigF(j,:) ) ;

    J(J==j) = [];
    C = [C,  1./ ( Z - Z(j) ) ];

    bigA = [];        JJ = [];
    %
    for ll = 1:ell
        mj = bigF(:,ll).*C - C*diag(smallF(:,ll));
        bigA = vertcat(bigA, mj);
        JJ = vertcat(JJ, M*(ll-1)+J(:) );
    end

    [~,~,V] = svd(bigA(JJ,:),0);   w = V(:,end);

    D  = C*w;
    bigR = bigF ;
    for ll = 1:ell
        mj = C*(w.*smallF(:,ll));
        bigR(J,ll) =  mj(J)./D(J);
    end

    %     check error
    err = norm(bigF(:) - bigR(:),inf);
    if err <  1e-14*norm(bigF(:),inf)
        fprintf('\n v_AAA achives the error %.4e and output # of support pts is %d \n',err,k)
        break,
    end

end

%  output index of support pts
mJ = setdiff((1:M).', J);




% % function handle
Rfun = cell(1,ell);

wj = w;    zj = z;
for ll = 1:ell
    fj1 = smallF(:,ll);
    R1 = @(zz) zz_ratevalB(zz, zj, fj1, wj);
    Rfun{ll} = R1;
end


% compute poles and zeros


zers = cell(1,ell);
for ll = 1:ell
    fj1 = smallF(:,ll);
    [pol,~, mj ] = zz_rootsB(zj,fj1, wj );
    zers{ll} = mj;
end




% ------------------------------------------------------------- 


if length(zj) < n + 1  
    n = length(zj) - 1;     
    fprintf('\n v_AAA compute degrees  nd = nk = %d   \n ', n  ) 
end









% % ---------------------------------------------------

%    lawson's iteration is used for intepolant form
% 
%         zz_v_aaa(t, bigF, n, 'lawson', maxit)
% 
%          solve:  
%                     |w1               |     |f1-r1|
%            min      |     w2          |     |f2-r2|
%                     |           :     |         :
%                     |               wk|     |fk-rk |
% 
%          by update :
%                    w_new = w_old*abs( error ) 
% 
%  in this code, support pts is deleted from the sample points 
% 
% % ---------------  lawson's iteration  -----------------


if nargin>3


    if varargin{1} == 'lawson'  &  isnumeric(varargin{2})

        maxit = varargin{2} + 1;       wA = w;
%         if maxit>1
%             fprintf('\n v_AAA with lawson iteration is used for v_aaa \n')
%         end


        bigS = 1;
        bigM0 = bigA(JJ,:);

        %         weigth for each fk
        wt = cell(1,ell);
        for ll = 1:ell
            wt{ll} = 1;
        end

        for kk = 1 : maxit
            %   bigS:  iteration matrix (weight matrix W)

            bigM = bigS.*bigM0;

            [~,~,V] = svd(bigM,0);    w = V(:,end);       % TLS by "economy SVD "
            %     [~,~,v] = svds(bigM,1,'smallest');
            %     v: minmium sigular vector

            % approximation of bigR ~= bigF
            D  = C*w;
            bigR = bigF ;
            for ll = 1:ell
                mj = C*(w.*smallF(:,ll));
                bigR(J,ll) =  mj(J)./D(J);
            end


            e1  = abs(bigF - bigR);


            %    ---------- check error and
            errA = norm(e1(:),"inf");
            if errA < 4e-15
                fprintf('\n v_AAA(%d) achives max error = %.4e and stop \n', kk-1, errA)
                wA = w;   kA = kk;
                break
            end

            %     remain the best "alp and bet"
            if errA < err
                wA = w;     kA = kk;
                err  = errA;
            end

            % compute error and update weight
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
            www = www(JJ);
            
            cst = norm(www,1);

            www = www(:)./cst;
            for l = 1:ell
                wt{l} = wt{l}/cst;
            end
            
            % sparse weight matrix
            %     bigS = sparse(1:ell*M,1:ell*M, sqrt( www(:) )  );
            bigS = sqrt(www(:));


        end

    %     --------------------------------------------------------------
 


%             ---------------  lawson's iteration  -----------------
%         output rational function and poles and zeros for lawson's result

        % % function handle
        Rfun = cell(1,ell);

        wj = wA;    zj = z;
        for ll = 1:ell
            fj1 = smallF(:,ll);
            R1 = @(zz) zz_ratevalB(zz, zj, fj1, wj);
            Rfun{ll} = R1;
        end


        % compute poles and zeros


        zers = cell(1,ell);
        for ll = 1:ell
            fj1 = smallF(:,ll);
            [pol,~, mj ] = zz_rootsB(zj,fj1, wj );
            zers{ll} = mj;
        end

 

    end


end






 