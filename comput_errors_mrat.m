function [RMSE, E, FnormErr] = comput_errors_mrat(Rfun,X,FX0) 




ell = size(FX0,2); 


mat_err = []; 


for k = 1:ell

    R1 = Rfun{k};   
    mat_err = [ mat_err,  R1(X) - FX0(:,k) ];

end



    npts = length(X); 
    RMSE  = norm(mat_err,"fro")/sqrt(npts); 
 

   E = max( sqrt( sum( mat_err.*conj(mat_err) , 2) ) ); 


   FnormErr = sqrt( sum( mat_err.*conj(mat_err) , 2) );  


