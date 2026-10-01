function check = test_block_aaa()

pts = linspace(0,10,20);
opts.tol = 1e-16;
opts.maxit = 10;
FF = @(z) [ 1/(z+1) 1/(z^2-5) z;
    1/(z^2+5+1i) (2+z^2)/(z^3 + 3*z^2 + 1) 7];
[R,rmse,out] = util_block_aaa(FF,pts,opts); 

check = [ min(rmse) norm(R(1) - FF(1)) ] < 1e-14;

end
