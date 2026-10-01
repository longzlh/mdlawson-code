%% The block AAA algorithm 
%  Ion Victor Gosea \and Stefan Guettel
%
%  July 2020
%
%  Tags: block AAA, BARYFUN

%% Introduction
% Since version 2.9 the RKToolbox provides two new utility functions 
% |util_block_aaa| and |util_rmse| for working with rational functions
% represented in a generalized barycentric form with matrix-valued weights, 
% 
% $$
%   R(z) = \left(\displaystyle \sum_{i=0}^d W_i/(z-z_i)\right)^{-1}
%          \left(\displaystyle \sum_{i=0}^d W_i F(z_i)/(z-z_i)\right).
% $$ 
%
% In this notation, the $W_i$ are $m\times m$ nonsingular weight matrices,
% the $z_i$ are pairwise distinct barycentric support points in the 
% complex plane, and the $F(z_i)$ are $m\times n$ matrices. 
% Different from the representation of rational functions computed by the
% RKFIT method [1], the AAA method [5] and its surrogate [2] and set-valued
% variants [4], the rational function $R(z)$ has a nonscalar "denominator 
% matrix polynomial." See [3] for more details and a number of numerical 
% experiments with different representations of rational matrix-valued
% functions.

%% Demonstration
% We focus on a simple $2\times 3$ matrix-valued function $F(z)$ to
% demonstrate the use of |util_block_aaa|. This function is defined as
% below and we sample it at 20 equidistant points on the interval $[0,10]$:

F = @(z) [ z/(z+1)      1/(z^2-5)  ;
             1/(z^2-5)    (2+z^2)/(z^3+3*z^2+1) ];
pts =    linspace(0, 11 ,  2000).';

%%
% Together with an error tolerance |tol| and a number of maximal iterations
% |maxit|, we have all that is needed to reapproximate $F$ using |util_block_aaa|:

opts.tol = 1e-12;
opts.maxit = 20;
[R,rmse,out] = util_block_aaa(F,pts,opts); 


% N = 100;
% A = gallery('tridiag', N);
% %A{2} = A{1};
% I = speye(N);
% F{1} = (A - 4*I)\(A - 2*I)^2/(A^2 - I);
% F{2} = (A - 4*I)\(A - 2*I)^3/(A^2 - I);
N = length(pts); 
Rs = zeros(4,N); 

for j = 1:N
    mj =  R(pts(j))   ; 
    Rs(:,j) =  mj(:)  ; 
end
Rs = Rs.'; 



ell = 4; 

FX = zeros(ell,N);   
for j = 1:N
    mj =  F(pts(j))   ; 
    FX(:,j) =  mj(:)  ; 
end
FX = FX.'; 


err_blkaaa  = norm(FX(:)-Rs(:),inf) 


% A = sparse(1:N,1:N,pts);  
 
for j = 1:ell
  Fmat{j} = util_build_real_matrix(FX(:, j));
end

A = util_build_real_matrix(pts); 

N = size(A,1); 

b = ones(N,1); 


xi = inf(1, 6);
% xi = inf(1,9);      % this gives the better results 
[xi,ratfun,misfit] = rkfit(Fmat, A, b, xi, 2, 1e-11 );

 for j = 1:ell
     rj = ratfun{j}; 
    mj =  rj(pts)   ; 
    R2(:,j) =  mj(:)  ; 
end
 
err_rkfit = norm(R2(:) - FX(:), inf)




%%
% The output |R| is a |baryfun| object and we can display its info as
% follows:

disp(R)

%%
% The output |rmse| stores the root mean squared error over all sampling
% points for each iteration of the block AAA method. In this case, we have
% resolved $F$ to about machine precision after 6 iterations:

semilogy(rmse)
axis tight, hold on
legend('block AAA')
xlabel('iteration'), ylabel('RMSE')
disp(['The final RMSE is ' num2str(util_rmse(pts,F,R)) ])

%% References
%
% [1] M. Berljafa and S. Guettel.
%     _The RKFIT algorithm for nonlinear rational approximation,_ 
%     SIAM J. Sci. Comput., 39(5):A2049--A2071, 2017.
%
% RKT_BIGBREAK
% 
% [2] S. Elsworth and S. Guettel. _Conversions between barycentric, RKFUN, 
%     and Newton representations of rational interpolants,_ 
%     Linear Algebra Appl., 576:246--257, 2019.
%
% RKT_BIGBREAK
%
% [3] I. V. Gosea and S. Guettel. _Algorithms for the rational 
%     approximation of matrix-valued functions,_ arXiv preprint 
%     2003.06410v1, 2020. (<https://arxiv.org/abs/2003.06410>)
%
% RKT_BIGBREAK
%
% [4] P. Lietaert, J. Perez, B. Vandereycken, and K. Meerbergen. _Automatic 
%     rational approximation and linearization of nonlinear eigenvalue 
%     problems,_ arXiv preprint 1801.08622, 2018.
%     (<https://arxiv.org/pdf/1801.08622.pdf>)
%
% RKT_BIGBREAK
%
% [5] Y. Nakatsukasa, O. Sete, and L. N. Trefethen. 
%     _The AAA algorithm for rational approximation,_
%     SIAM J. Sci. Comput., 40(3):A1494--A1522, 2018.
