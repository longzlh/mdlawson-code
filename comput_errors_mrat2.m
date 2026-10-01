function [RMSE, E, FnormErr] = comput_errors_mrat2(Rfun,X,FF0,ell) 

npts = length(X);   

% ell = 4 ;  

 

%     err1 = zeros(1,npts);
mat_err =  zeros(npts, ell);
    % compute max F-error
    for k = 1:npts
        mj = FF0{k};
        mj2 = Rfun(X(k) );
        mat_err(k,:) = ( mj(:) - mj2(:) ).';
    end

% figure; plot(X./1i, err1,'LineWidth',2); grid on,
    RMSE = norm(mat_err,"fro"); 
    E = max( sqrt( sum( mat_err.*conj(mat_err) , 2) ) );


       FnormErr = sqrt( sum( mat_err.*conj(mat_err) , 2) );  