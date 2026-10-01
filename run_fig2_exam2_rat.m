clear,
close all,
%
% Example 5.5 from Ion Victor Gosea & Stefan Guttel
% Algorithms for the rational approximation of matrix-valued functions

% % test vector value function

% F = @(z) [z.*(1-2*z.*cot(2*z))./(tan(z) - z) + 10, ...
%     z.*(2*z - sin(2*z))./ ( sin(2*z).*tan(z) - z ), ...
%     z.*(2*z - sin(2*z))./ ( sin(2*z).*tan(z) - z ), ...
%     z.*(1-2*z.*cot(2*z))./(tan(z) - z) + 4 ];


F = @(z) [z.*(1-2*z.*cot(2*z))./(tan(z) - z) + 10, ...
    z.*(2*z - sin(2*z))./ ( sin(2*z).*(tan(z) - z) ), ...
    z.*(2*z - sin(2*z))./ ( sin(2*z).*(tan(z) - z) ), ...
    z.*(1-2*z.*cot(2*z))./(tan(z) - z) + 4 ];



% sample pts
M =  500 ;
X = logspace(-2, 1, M).'*1i;
% X = linspace(0.01,10,M).'*1i;

ell = 4;      % sample date  FX = [f1(X), f2(X),..., f_ell(X)]
FX = X*zeros(1,ell);
for j = 1:M
    xj = X(j);    mj = F(xj);
    FX(j,:) = mj(:).';
end




%      --------------- vec fun -------------------------
%
%     rational approximation :
%               min | f1-p1/q | + |f2-p2/q| + ...+ |fk - pk/q|
%
nd = 6;               % degree of q

% deg_vec = zeros(1,ell) + 8;        % (nd = nk)

deg_vec =  -1 + [7 7 7 7];                 % (nd~=nk)

maxit = 20;


tic

% ------  aaa method for vector-valued function

% [Rfun, mJ, J,  pol,  zers ] = zz_v_aaa_lawson(X, FX, nd, deg_vec, maxit);
% %
% err = []; dk = [];


% ---- d-lawson method for vector-valued function

% [Rfun, relgap, dk, err, pol, zers, w] ...
%     = zz_v_dlawson(X, FX, nd, deg_vec, maxit );

[Rfun, relgap, dk, err, pol, zers, w] ...
    = m_d_lawson(X, FX, nd, deg_vec, maxit );

% [RMSE, ER] = comput_errors_mrat(Rfun,X, FX) ;
% [Rfun, Ak, relgap, dk, err, zers ] ...
%     = zz_v_bestpoly(X, FX, deg_vec, maxit );


toc




R1 = Rfun{1};
R2 = Rfun{2};
R3 = Rfun{3};
R4 = Rfun{4};

err1 = norm(R1(X) - FX(:,1), "inf");
err2 = norm(R2(X) - FX(:,2), "inf");
err3 = norm(R3(X) - FX(:,3), "inf");
err4 = norm(R4(X) - FX(:,4), "inf");

E = [err1,err2,err3,err4];
fprintf('\n Max Errors of  r1,r2,r3 , r4 \n  %e, %e, %e, %e \n',E)


ERRS = [];
for k = 1:ell
    r = Rfun{k};
    mj = abs( r(X) - FX(:,k) );
    ERRS = [ERRS, mj ];
end
RMSE = norm(ERRS,'fro')/sqrt(M),

% RMSE = 0 ;
% for k = 1:ell
%     RX = Rfun{ell};
%     mj =  ( RX(X) - FX(:,ell) );
%     RMSE = RMSE + mj'*mj ;
% end
% RMSE =  sqrt(RMSE/M) ,

% ----------------- plot error curves ----------------------


fig1 = figure ;       LW = 'LineWidth';
subplot(2,3,1),
res = ERRS(:,1);
plot(X./1i, abs(res),'-',LW,2); grid on,
xlabel('$ x_{\ell}\in [10^{-2},10]{\rm i} $','Interpreter','latex');
ylabel('$|f_{11}(x_{\ell}) - r_{11}(x_{\ell})|$','Interpreter','latex')
str = { ' Error curve of r_{11}' } ;   title(str)
set(gca,'FontSize',14)
% ylim([0,4e-3])

subplot(2,3,2),
res = ERRS(:,2);
plot(X./1i, abs(res),'-',LW,2); grid on,
xlabel('$ x_{\ell}\in [10^{-2},10]{\rm i} $','Interpreter','latex');
ylabel('$|f_{12}(x_{\ell}) - r_{12}(x_{\ell})|$','Interpreter','latex')
str = { ' Error curve of r_{12}' } ;   title(str)
set(gca,'FontSize',14)
% ylim([0,1.8e-3])


subplot(2,3,4),
res = ERRS(:,3);
plot(X./1i, abs(res),'-',LW,2); grid on,
xlabel('$ x_{\ell}\in [10^{-2},10]{\rm i} $','Interpreter','latex');
ylabel('$|f_{21}(x_{\ell}) - r_{21}(x_{\ell})|$','Interpreter','latex')
str = { ' Error curve of r_{21}' } ;   title(str)
set(gca,'FontSize',14)
% ylim([0,1.8e-3])

subplot(2,3,5),
res = ERRS(:,4);
plot(X./1i, abs(res),'-',LW,2); grid on,
xlabel('$ x_{\ell}\in [10^{-2},10]{\rm i} $','Interpreter','latex');
ylabel('$ |f(x_{\ell})_{22} - r_{22}(x_{\ell})|$','Interpreter','latex')
str = { ' Error curve of r_{22}' } ;   title(str)
set(gca,'FontSize',14)
% ylim([0,1.8e-3])



subplot(2,3,6),
plot( err,'-',LW,2); hold on,
plot(dk,'-',LW,2);   grid on,
ylim(err(end)*[0,3] )
s = legend(' $ \sqrt{e(R^{(k)}})$','$ \sqrt{d(\mathbf{w}^{(k)})}$');
set(s,'Interpreter','latex');
xlabel('Lawson iteration k');
set(gca,'FontSize',14)






subplot(2,3,3)
plot(X./1i,  sqrt( sum( ERRS.^2, 2) ) , 'LineWidth', 2) ; grid on ,
hold on,
maxF = norm( sqrt( sum( ERRS.^2, 2) ),  'inf' );
plot(X./1i, maxF + 0*X, 'r-', 'LineWidth', 2);
complement = max(abs(w.*(maxF.^2-( sum( ERRS.^2, 2)))))
xlabel('$ x_{\ell}\in [10^{-2},10]{\rm i} $','Interpreter','latex');
ylabel('$ \| F(x_{\ell}) - R(x_{\ell})\|_F $ ','Interpreter','latex');
set(gca,'FontSize',14)
ylim([0,maxF*1.1])

%str = {'complementary =',complement};
str = {'type (6, 6)'};
text(2, 0.8e-6,  str, 'fontsize',22,'Color','r','Interpreter','latex');
set(fig1,'InnerPosition', [600 600 1200 500])


saveas(fig1,'exam2_rat66_errs','epsc')



