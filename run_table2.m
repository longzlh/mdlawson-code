function [E,RMSE,T] = run_table2()
%RUN_TABLE2  Example 8.2 data: max error, RMSE, and runtime.
% Example 5.1 from Ion Victor Gosea & Stefan Guttel
% Algorithms for the rational approximation of matrix-valued functions
% pts = 2*exp(2i*pi*(1:npts)/npts)';        % sampling points


npts =  500;
pts = logspace(-2,1,npts)'*1i;


F = @(z) [z.*(1-2*z.*cot(2*z))./(tan(z) - z) + 10, ...
    z.*(2*z - sin(2*z))./ ( sin(2*z).*(tan(z) - z) ), ...
    z.*(2*z - sin(2*z))./ ( sin(2*z).*(tan(z) - z) ), ...
    z.*(1-2*z.*cot(2*z))./(tan(z) - z) + 4 ];


ell = 4;

noiselevel = [0 1e-10 1e-8 1e-6 1e-4];
E_md =[];   E_va =[];     E_ba =[];
E_lw =[];   E_rk = [];    E_vf =[];
RMSE_md =[];   RMSE_va =[];     RMSE_ba =[];
RMSE_lw =[];   RMSE_rk = [];    RMSE_vf =[];
T_md = [];  T_va = [];  T_ba = [];
T_lw = [];  T_rk = [];  T_vf = [];
%noise on the matrix-valued fun

X = pts;

deg = 10;
nd = deg;

deg_vec = zeros(1,ell) + nd ;

maxit = 10 ;


for mm = 1:length(noiselevel)
    sig = noiselevel(mm);
    %% prepare data
    %   the length of function F

    % % sample date  FX = [f1(X), f2(X),..., f_ell(X)]

    FX = zeros(ell, npts );
    for j = 1:npts
        xj = X(j);       mj = F(xj);
        FX(:,j) = mj(:);
    end

    FX0 = FX.';

    FX = FX0 + sig* ( randn(size(FX0 ) ) + 1i*randn(size(FX0) ) ) ;


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



    tic_md = tic;

    %Rfun = zz_v_dlawson(X, FX, nd, deg_vec, maxit);
    Rfun = m_d_lawson(X, FX, nd, deg_vec, maxit);
    T_md = [T_md; toc(tic_md)];

    [rmse_md, e_md] = comput_errors_mrat(Rfun,X, FX0) ;
    E_md = [E_md;  e_md];
    RMSE_md  = [RMSE_md;  rmse_md];

    %% Vector-AAA
    %     deg_vec = [nn, nt, nt, nr] ;


    tic_va = tic;
    Rfun = zz_v_aaa_lawson(X, FX, nd, deg_vec, maxit);
    T_va = [T_va; toc(tic_va)];

    [rmse_va,  e_va] = comput_errors_mrat(Rfun,X,FX0) ;

    E_va = [E_va;  e_va];
    RMSE_va = [RMSE_va;  rmse_va];

    %%  block aaa,

    opts.tol = 0;
    opts.maxit = nd ;
    opts.return = 'last';
    tic_ba = tic;
    R_ba = block_aaa(FF,pts,opts);
    T_ba = [T_ba; toc(tic_ba)];

    [rmse_ba, e_ba] = comput_errors_mrat2(R_ba,X,FF0,ell) ;
    E_ba = [E_ba;  e_ba];
    RMSE_ba = [RMSE_ba;  rmse_ba];


    %%   RKFIT
    ops.iter = maxit;
    ops.tol = 0;

    deg = nd;

    xi = inf(deg,1);

    tic_rkfit = tic;
    [Rkfit,rmse_rkfit] = baryfit6(FF,pts,xi,ops);
    % xi -- initial poles
    T_rk = [T_rk; toc(tic_rkfit)];


    [rmse_rk, e_rk] = comput_errors_mrat2(Rkfit,X,FF0,ell) ;
    E_rk = [E_rk;  e_rk];
    RMSE_rk = [RMSE_rk;  rmse_rk];

    %% Loewner


    deg = nd;

    tic_lw = tic;
    R_loew  = loewfit3(FF,pts,deg);
    T_lw = [T_lw; toc(tic_lw)];
    % xi -- initial poles


    [rmse_lw, e_lw] = comput_errors_mrat2(R_loew,X,FF0,ell) ;

    E_lw = [E_lw;  e_lw];
    RMSE_lw = [RMSE_lw;  rmse_lw];

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

    tic_vf = tic;
    R_vf = baryfit5(FF,pts,deg,opts);
    T_vf = [T_vf; toc(tic_vf)];



    [rmse_vf, e_vf] = comput_errors_mrat2(R_vf,X,FF0,ell) ;

    E_vf = [E_vf;  e_vf];
    RMSE_vf = [RMSE_vf;  rmse_vf];


end


RMSE = [RMSE_md, RMSE_va, RMSE_ba,  RMSE_rk, RMSE_lw, RMSE_vf];

E = [E_md, E_va, E_ba,  E_rk, E_lw, E_vf];
T = [T_md, T_va, T_ba, T_rk, T_lw, T_vf];

% comput_errors_mrat already returns RMSE = ||err||_F / sqrt(m).
% comput_errors_mrat2 returns ||err||_F. Divide by sqrt(m) so every
% column matches (eq:RMSE) used in Table 8.2.
RMSE(:, 3:6) = RMSE(:, 3:6) / sqrt(npts);

end
