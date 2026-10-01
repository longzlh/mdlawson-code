clear,
close all,



npts =  500;
pts = logspace(-2,1,npts)'*1i;


F = @(z) [z.*(1-2*z.*cot(2*z))./(tan(z) - z) + 10, ...
        z.*(2*z - sin(2*z))./ ( sin(2*z).*(tan(z) - z) ), ...
        z.*(2*z - sin(2*z))./ ( sin(2*z).*(tan(z) - z) ), ...
        z.*(1-2*z.*cot(2*z))./(tan(z) - z) + 4 ];


ell = 4; 

 
 
%noise on the matrix-valued fun

X = pts;

deg = 6;
nd = deg;

deg_vec = zeros(1,ell) + nd ;

maxit = 20 ;


    %% prepare data
    %   the length of function F
 
    % % sample date  FX = [f1(X), f2(X),..., f_ell(X)]

    FX = zeros(ell, npts );
    for j = 1:npts
        xj = X(j);       mj = F(xj);
        FX(:,j) = mj(:);
    end

    FX0 = FX.';

    FX = FX0 ;


    FF = cell(1,npts);     FF0 = FF;
    for j = 1:npts
        mj = FX(j,:);
        FF{j} = reshape(mj,2,2);

        mj = FX0(j,:);
        FF0{j} = reshape(mj,2,2);

    end

    %  FX   perturb,   FX0  original,
    FX0 = FX;   FF0 = FF;

    %%     m-d-lawson
 fig1 = figure;  LW = 'LineWidth';     

    %Rfun = zz_v_dlawson(X, FX, nd, deg_vec, maxit);
    Rfun = m_d_lawson(X, FX, nd, deg_vec, maxit);

    [~, ~, FnormErr_md] = comput_errors_mrat(Rfun,X, FX0) ;
 
subplot(2,3,1);  plot(X./1i,FnormErr_md,LW,2); grid on, 
hold on,  maxF = norm( FnormErr_md,  'inf' ); 
plot(X./1i, maxF + 0*X, 'r-', 'LineWidth', 2)  
ylabel('$ \| F(x_{\ell}) - R(x_{\ell}) \|_{F} $','Interpreter','latex'); 
xlabel('$ x_{\ell}\in [10^{-2},10]{\rm i} $','Interpreter','latex');  
set(gca,'FontSize',14); 
title('m-d-Lawson(20)')
str = {'type (6, 6)'};
text(2, 1e-6,  str, 'fontsize',22,'Color','r','Interpreter','latex');
    %% Vector-AAA
%     deg_vec = [nn, nt, nt, nr] ;


    Rfun = zz_v_aaa_lawson(X, FX, nd, deg_vec, maxit);

    [~,~,  FnormErr_va] = comput_errors_mrat(Rfun,X,FX0) ;
 
subplot(2,3,2);  plot(X./1i,FnormErr_va,LW,2); grid on, 
hold on,  maxF = norm( FnormErr_va,  'inf' ); 
plot(X./1i, maxF + 0*X, 'r-', 'LineWidth', 2)  
ylabel('$ \| F(x_{\ell}) - R(x_{\ell}) \|_{F} $','Interpreter','latex'); 
xlabel('$ x_{\ell}\in [10^{-2},10]{\rm i} $','Interpreter','latex');  
set(gca,'FontSize',14); 
title('v-AAA-Lawson(20)')
str = {'type (6, 6)'};
text(2, 2e-6,  str, 'fontsize',22,'Color','r','Interpreter','latex');
    %%  block aaa,

    opts.tol = 0;
    opts.maxit = deg;
    opts.return = 'last';
    tic_ba = tic;

    R_ba = block_aaa(FF,pts,opts);
    toc(tic_ba);

     [~,~,  FnormErr_ba] = comput_errors_mrat2(R_ba,X,FF0,ell) ;
 

subplot(2,3,3);  plot(X./1i,FnormErr_ba,LW,2); grid on, 
hold on,  maxF = norm( FnormErr_ba,  'inf' ); 
plot(X./1i, maxF + 0*X, 'r-', 'LineWidth', 2)  
ylabel('$ \| F(x_{\ell}) - R(x_{\ell}) \|_{F} $','Interpreter','latex'); 
xlabel('$ x_{\ell}\in [10^{-2},10]{\rm i} $','Interpreter','latex');  
set(gca,'FontSize',14); 
title('block-AAA')
str = {'type (6, 6)'};
text(2, 2e-6,  str, 'fontsize',22,'Color','r','Interpreter','latex');
    %%   RKFIT
    ops.iter = maxit;
    ops.tol = 0;

    deg = nd;

    xi = inf(deg,1);

    tic_rkfit = tic;
    [Rkfit,rmse_rkfit] = baryfit6(FF,pts,xi,ops);
    % xi -- initial poles
    toc(tic_rkfit)


    [~,~,  FnormErr_rk] = comput_errors_mrat2(Rkfit,X,FF0,ell) ;
 
subplot(2,3,4);  plot(X./1i,FnormErr_rk,LW,2); grid on, 
hold on,  maxF = norm( FnormErr_rk,  'inf' ); 
plot(X./1i, maxF + 0*X, 'r-', 'LineWidth', 2)  
ylabel('$ \| F(x_{\ell}) - R(x_{\ell}) \|_{F} $','Interpreter','latex'); 
xlabel('$ x_{\ell}\in [10^{-2},10]{\rm i} $','Interpreter','latex');  
set(gca,'FontSize',14); 
title('RKFIT(20)')
str = {'type (6, 6)'};
text(2, 5e-6,  str, 'fontsize',22,'Color','r','Interpreter','latex');
    %% Loewner


    deg = nd;

    R_loew  = loewfit3(FF,pts,deg);
    % xi -- initial poles


 [~,~,  FnormErr_lw] = comput_errors_mrat2(R_loew,X,FF0,ell) ;

 
subplot(2,3,5);  plot(X./1i,FnormErr_lw,LW,2); grid on, 
hold on,  maxF = norm( FnormErr_lw,  'inf' ); 
plot(X./1i, maxF + 0*X, 'r-', 'LineWidth', 2)  
ylabel('$ \| F(x_{\ell}) - R(x_{\ell}) \|_{F} $','Interpreter','latex'); 
xlabel('$ x_{\ell}\in [10^{-2},10]{\rm i} $','Interpreter','latex');  
set(gca,'FontSize',14); 
title('Loewner')
str = {'type (6, 6)'};
text(2, 7e-2,  str, 'fontsize',22,'Color','r','Interpreter','latex');
    %%  Vector Fitting
    opts.Niter2 = 5;
    opts.cmplx_ss = 1;
    opts.poletype='linlogcmplx';
    opts.Niter1=0;
    opts.asymp=2;
    opts.plot=0;
    opts.cmplx_ss=1;
    opts.stable = 0;
    opts.weightparam=1;
    opts.screen=0;

    deg = nd;

    R_vf = baryfit5(FF,pts,deg,opts);



 [~,~,  FnormErr_vf] = comput_errors_mrat2(R_vf,X,FF0,ell) ;

 
subplot(2,3,6);  plot(X./1i,FnormErr_vf,LW,2); grid on, 
hold on,  maxF = norm( FnormErr_vf,  'inf' ); 
plot(X./1i, maxF + 0*X, 'r-', 'LineWidth', 2)  
 ylabel('$ \| F(x_{\ell}) - R(x_{\ell}) \|_{F} $','Interpreter','latex'); 
xlabel('$ x_{\ell}\in [10^{-2},10]{\rm i} $','Interpreter','latex');  
 set(gca,'FontSize',14); 
title('VF')
str = {'type (6, 6)'};
text(2, 1e-3,  str, 'fontsize',22,'Color','r','Interpreter','latex');
% % output figs

set(fig1,'InnerPosition', [600 600 1200 500])

saveas(fig1,'exam2_rat66_compares','epsc')




