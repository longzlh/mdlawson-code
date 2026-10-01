 clear,
close all

% % % % % % % % % % % % % %
 
data_ieee,
 
%% random perturb 
    sig = 1e-8;
 
    mj1 = randn(size(y1) ) + 1i*randn(size(y1) )   ;     y1 = y1 + sig*mj1;
    mj2 = randn(size(y2) ) + 1i*randn(size(y2) )   ;     y2 = y2 + sig*mj2;
    mj3 = randn(size(y3) ) + 1i*randn(size(y3) )   ;     y3 = y3 + sig*mj3;
    % #0 è®°å???å§??°æ??
    t0 = t;  y10 = y1;  y20 = y2;  y30 = y3;
 
 
    %% Our Algorithms
    pts = t;   Ns = length(pts);

    deg_vec = [nn, nt, nr];

    maxit = 10;

    X = pts;      FX = [y1,y2,y3];


    [Rfun, mJ, J, pol, zers, err_va] = zz_v_aaa_lawson(X, FX, nd, deg_vec, maxit);

    %[Rfun2 , relgap, dk, err, pol2, zers2 ] = zz_v_dlawson(X, FX, nd, deg_vec, maxit);
    [Rfun2 , relgap, dk, err, pol2, zers2 ] = m_d_lawson(X, FX, nd, deg_vec, maxit);



% %     vec-aaa
    R1 = Rfun{1}; R2 = Rfun{2}; R3 = Rfun{3};

    zer_D1 = pol(1:nd);
    zer_N1 = zers{1};  zer_T1 = zers{2};  zer_R1 = zers{3};
 
    % % sort the zeros
     zer_D1 = comput_2vecs_err(zer_D,zer_D1); 
     zer_N1 = comput_2vecs_err(zer_N,zer_N1); 
     zer_R1 = comput_2vecs_err(zer_R,zer_R1); 
     zer_T1 = comput_2vecs_err(zer_T,zer_T1); 
    %
    errN = norm(zer_N1-zer_N)/norm(zer_N);
    errT = norm(zer_T1-zer_T)/norm(zer_T);
    errR = norm(zer_R1-zer_R)/norm(zer_R);
    errD = norm(zer_D1-zer_D)/norm(zer_D);

 
    relE = [errD, errN, errT, errR];
 
% % m-d-lawson 
    mR1 = Rfun2{1}; mR2 = Rfun2{2}; mR3 = Rfun2{3};

   zer_D2 = pol2(1:nd);
   zer_N2 = zers2{1};  
   zer_T2 = zers2{2};  
   zer_R2 = zers2{3};
 
    % % sort the zeros
     zer_D2 = comput_2vecs_err(zer_D,zer_D2); 
     zer_N2 = comput_2vecs_err(zer_N,zer_N2); 
     zer_R2 = comput_2vecs_err(zer_R,zer_R2); 
     zer_T2 = comput_2vecs_err(zer_T,zer_T2); 
    %
    errN2 = norm(zer_N2-zer_N)/norm(zer_N);
    errT2 = norm(zer_T2-zer_T)/norm(zer_T);
    errR2 = norm(zer_R2-zer_R)/norm(zer_R);
    errD2 = norm(zer_D2-zer_D)/norm(zer_D);

 
    relE2 = [errD2, errN2, errT2, errR2];


%% output table
str = { 'error of roots p_D', 'error of roots p_N',  'error of roots p_T','error of roots p_R'};
 EE = [relE; relE2],
% array2table(EE ,  'VariableNames', str ),



LW = 'LineWidth';       FS = 'FontSize';    MK = 'MarkerSize'; 
fig1 = figure(1);
subplot(2,3,1)
plot(zer_D1,'bo',MK,9);hold on
plot(zer_D, 'r.',MK,22);
title('original and computed zeros of p_D',FS,14); grid on,  
legend('computed','original')

subplot(2,3,2)
plot(zer_N1,'bo',MK,9);hold on
plot(zer_N,'r.', MK,22);
title('original and computed zeros of p_N',FS,14); grid on,  
legend('computed','original')

subplot(2,3,4)
plot(zer_T1,'bo',MK,9);hold on
plot(zer_T,'r.', MK,22);
title('original and computed zeros of p_T',FS,14); grid on,   
legend('computed','original')

subplot(2,3,5)
plot(zer_R1,'bo',MK,9);hold on
plot(zer_R,'r.', MK,22);
title('original and computed zeros of p_R',FS,14); grid on,   
legend('computed','original')
% set(gcf,"Units",'normalized','InnerPosition',[.03,.1,.6,.75])
% set(gca,FS,16)

% plot the fitting curve
subplot(2,3,[3,6])
t_new = t0;
plot_zlh3(t0,y10,y20,y30,t_new,R1(t_new),R2(t_new),R3(t_new));
% set(gcf,"Units",'normalized','InnerPosition',[.03,.1,.5,.6])
set(gca,FS,16)

str = {' v-AAA-Lawson(10)   ','$\zeta = 10^{-8}$'}; 
text(-1.6, -80,  str, 'fontsize',18,'Color','r','Interpreter','latex' );

%---------------------------  


fig2 = figure(3);
subplot(2,3,1)
plot(zer_D2,'bo',MK,9);hold on
plot(zer_D,'r.', MK,22);
title('original and computed zeros of p_D',FS,14); grid on,  
legend('computed','original')

subplot(2,3,2)
plot(zer_N2,'bo',MK,9);hold on
plot(zer_N,'r.', MK,22);
title('original and computed zeros of p_N',FS,14); grid on,  
legend('computed','original')

subplot(2,3,4)
plot(zer_T2,'bo',MK,9);hold on
plot(zer_T,'r.', MK,22);
title('original and computed zeros of p_T',FS,14); grid on,  
legend('computed','original')

subplot(2,3,5)
plot(zer_R2,'bo',MK,9);hold on
plot(zer_R,'r.', MK,22);
title('original and computed zeros of p_R',FS,14); grid on,   
legend('computed','original')
% set(gcf,"Units",'normalized','InnerPosition',[.03,.1,.6,.75])
% set(gca,FS,16)

% plot the fitting curve
subplot(2,3,[3,6])
t_new = t0;
plot_zlh3(t0,y10,y20,y30,t_new,mR1(t_new),mR2(t_new),mR3(t_new));
% set(gcf,"Units",'normalized','InnerPosition',[.03,.1,.5,.6])
set(gca,FS,16)
str = {' m-d-Lawson(10)   ','$\zeta = 10^{-8}$'}; 
text(-1.6, -80,  str, 'fontsize',18,'Color','r','Interpreter','latex' );


set(fig1,'InnerPosition', [600 600 1200 400])
set(fig2,'InnerPosition', [1000 600 1200 400])
% saveas 

% fig_mdlawson_HW_nosie8

saveas(fig2, 'fig_mdlawson_HW_nosie8','epsc')
saveas(fig1, 'fig_vAAAlawson_HW_noise8','epsc')
% % ===================== plot_figure ===========================

function fig0 = plot_zlh3(t0,y10,y20,y30, tr,R1,R2,R3)  

   if ~isreal(t0) 
      t0 = imag(t0);  tr = imag(tr);  
   end
 
fig0 =  plot(t0, 20*log10(abs(y10) ),'b-','LineWidth',1.5); hold on, 
  plot(t0, 20*log10(abs(y20) ),'r-','LineWidth',1.5);   
  plot(t0, 20*log10(abs(y30) ),'m-','LineWidth',1.5);   
%   
 plot(tr, 20*log10(abs(R1)),'k-.','LineWidth',1);  
 plot(tr, 20*log10(abs(R2)),'k-.','LineWidth',1); 
  plot(tr, 20*log10(abs(R3)),'k-.','LineWidth',1); hold off,
 grid on, 
 s = legend('$ y_{1}$', '$ y_{2}$ ','$ y_{3}$ ' );
 set(s,'Interpreter','latex','location','best')
 grid on, xlabel('Normalized Frequency'); 
%  title('ND + 1 - 1,\; TR + 0 ')
 ylabel('dB');    axis([t0(1) , t0(end)  , -180,20]);
% end

end






% ----------------- subroutine -------------------------


function [x1,err] = comput_2vecs_err(x,y) 

lx = length(x); 
err = zeros(lx,1); 
x1 = err; 

for k = 1:lx
       xk = x(k) ;
       [eI, I] = min(abs(y - xk)); 
       err(k) = eI; 
       x1(k) = y(I); 
end

end
