clear,    close all,
%
% Example 5.5 from Ion Victor Gosea & Stefan Guttel
% Algorithms for the rational approximation of matrix-valued functions

% % test vector value function
% lam = 1;
F = @(z) [z.*(1-2*z.*cot(2*z))./(tan(z) - z) + 10, ...
    z.*(2*z - sin(2*z))./ ( sin(2*z).*(tan(z) - z) ), ...
    z.*(2*z - sin(2*z))./ ( sin(2*z).*(tan(z) - z) ), ...
    z.*(1-2*z.*cot(2*z))./(tan(z) - z) + 4 ];

% sample pts
M = 500;
X = logspace(-2, 1, M).'*1i;

% M = 4000;
% X = linspace(0.01,10, M).'*1i;


% the length of function F
ell = 4;
% sample date  FX = [f1(X), f2(X),..., f_ell(X)]
FX = X*zeros(1,ell);
for j = 1:M
    xj = X(j);    mj = F(xj);
    FX(j,:) = mj(:).';
end

% sig = 1e-4;
% FX = FX + sig*rand(size(FX));

%% vec_aaa  and m_dlawson


% deg_vec = [5,  5,  5   ] ;
deg_vec0  = zeros(1,ell) ;

% if length(deg_vec)~= ell
%     error('given data and the numer of degrees not match!');
% end

maxit = 30;


RMSE_va = [];   RMSE_md = [];
E_va = [];        E_md = [];


RMSE2_va = [];   RMSE2_md = [];
E2_va = [];        E2_md = [];

nd = 20;    nk = nd;

for ijk = 1:nk

    deg_vec = deg_vec0 + ijk;

    Rfun = zz_v_aaa_lawson(X, FX,   nd,  deg_vec, maxit);
    

    %Rfun2 = zz_v_dlawson(X,    FX,    nd,    deg_vec,  maxit);
    Rfun2 = m_d_lawson(X, FX,   nd,  deg_vec, maxit);

    [rmse_va, e_va] = comput_errors_mrat(Rfun,X, FX) ;

    [rmse_md, e_md] = comput_errors_mrat(Rfun2,X, FX) ;


    RMSE_va = [RMSE_va; rmse_va];
    RMSE_md = [RMSE_md; rmse_md];

    E_va = [E_va; e_va];
    E_md = [E_md; e_md];




    % % ======================

    nd2 = ijk + 2;

    Rfun3 = zz_v_aaa_lawson(X, FX, nd2,    deg_vec, maxit);

    Rfun4 = zz_v_dlawson(X,    FX,    nd2,    deg_vec,  maxit);


    [rmse2_va, e2_va] = comput_errors_mrat(Rfun3,X, FX) ;

    [rmse2_md, e2_md] = comput_errors_mrat(Rfun4,X, FX) ;


    RMSE2_va = [RMSE2_va; rmse2_va];
    RMSE2_md = [RMSE2_md; rmse2_md];

    E2_va = [E2_va; e2_va];
    E2_md = [E2_md; e2_md];



end



fig3 = figure ;

subplot(2,2,1);
semilogy(1:nk,RMSE_va, 'r.-','MarkerSize', 26);
hold on,
semilogy(1:nk,RMSE_md, 'bo-','MarkerSize', 9);
grid on
xlabel('$ n_{ij} $ ','Interpreter','latex','FontSize',20)
ylabel('RMSE')
str = {'$ d = 20 $ '};
text(2, 1e-10,  str, 'fontsize',22,'Color','b','Interpreter','latex');
title('RMSE of v-AAA-Lawson(20) and m-d-Lawson(20) ');
ylim([1e-15,10])
legend('v-AAA-Lawson','m-d-Lawson')



subplot(2,2,2);
semilogy(1:nk,E_va, 'r.-','MarkerSize', 26);
hold on,
semilogy(1:nk,E_md, 'bo-','MarkerSize', 9);
grid on
xlabel('$ n_{ij}$ ','Interpreter','latex','FontSize',20)
ylabel('$\sqrt{e(R)}$','Interpreter','latex','FontSize',12)
str = {'$ d = 20 $ '};
text(2, 1e-10,  str, 'fontsize',22,'Color','b','Interpreter','latex');
title(' Max error of v-AAA-Lawson(20) and m-d-Lawson(20) ');
ylim([1e-15,100])
legend('v-AAA-Lawson','m-d-Lawson')



subplot(2,2,3);
semilogy(1:nk,RMSE2_va, 'r.-','MarkerSize', 26);
hold on,
semilogy(1:nk,RMSE2_md, 'bo-','MarkerSize', 9);
grid on
xlabel('$ n_{ij}$ ','Interpreter','latex','FontSize',20)
ylabel('RMSE')
str = {'$ d = n_{ij} + 2$ '};
text(2, 1e-10,  str, 'fontsize',22,'Color','b','Interpreter','latex');
title('RMSE of v-AAA-Lawson(20) and m-d-Lawson(20) ');
ylim([1e-15,10])
legend('v-AAA-Lawson','m-d-Lawson')



subplot(2,2,4);
semilogy(1:nk,E2_va, 'r.-','MarkerSize', 26);
hold on,
semilogy(1:nk,E2_md, 'bo-','MarkerSize', 9);
grid on
xlabel('$ n_{ij}$ ','Interpreter','latex','FontSize',20)
ylabel('$\sqrt{e(R)}$','Interpreter','latex','FontSize',12)
str = {'$ d = n_{ij} + 2$ '};
text(2, 1e-10,  str, 'fontsize',22,'Color','b','Interpreter','latex');
title(' Max error of v-AAA-Lawson(20) and m-d-Lawson(20) ');
ylim([1e-15,10])
legend('v-AAA-Lawson','m-d-Lawson')

set(fig3,'InnerPosition', [1000 600 1000 600])

saveas(fig3, 'ex2_lawson20_converg','epsc')
