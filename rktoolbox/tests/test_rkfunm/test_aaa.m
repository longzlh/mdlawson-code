function check = test_aaa()

fun = @(z) expm(z*[-1,1;.5 -1]);
Z = linspace(-4-40i, -4+40i, 100);
[rat,pol,res,zer,z,f,w,errvec] = util_aaa(fun,Z,1e-12);
ratm = util_aaa(fun,Z,1e-12);
for j = 1:length(Z),
   err(j) = norm(fun(Z(j)) - ratm(Z(j)))/norm(fun(Z(j))); 
end
check = max(err) < 1e-8;

end
