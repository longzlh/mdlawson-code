function [x, An] = zz_roots(H, d)
%
% comput roots of p = PHI*d, where PHI = [phi0,...,phi_n]
% is Vander type basis matrix
% from Arnoldi process with hessenberg H
% input:  hessenberg matrix H
%              coeffs d that p(z) = polyvalA(d,H,z)
%
% output: x // zeros of p(x)
%              An // leading coeff of p(x)
 

II = find(abs(d)>1e-15,1,'last');
if isempty(II)
    error('Input coeffs dk are all zeros');
elseif II>1
    d = d(1:II);
    H = H(1:II,1:II-1);
end


if  II <=1
    x = [];   An = []; 
    return
end


[m,n] = size(H);
if m ~= n+1, error('Input Hessberg matrix is wrong!'); end


An = d(end)./prod( diag(H,-1) );


% meth = 1;       % qr method
meth = 2;         % qz for general eigsolve



d = d./norm(d,2);

switch meth
    case 1

        %         may ill-conditon if d(n+1) is small
        d1 = d./d(n+1) ;

        H1 = H(1:n,1:n);


        cc = d1(1:n).*H(m,n);

        H1(:,n) = H1(:,n) - cc;

        x = eig(H1);

        %         x = sort(x);
        %         [~,x] = eig(H1);
        %         x = sort(diag(x));

    case 2
        % general eig : //   S^-1*A = H1
        A = H(1:n,1:n);
        S = eye(n);         S(n,n) = d(n+1);

        A(:,n) =  A(:,n)*d(n+1) - H(m,n)*d(1:n);

        [~,x] = eig(A, S, 'qz' );
        x = sort(diag(x));
        %         x = eig(A, S, 'qz' );
end

x = sort(x);


% x = x( ~isinf(x) ) ;



opts.newton = 0;
% ---- newton method for more accurate solution
% ---- one more newton step
if opts.newton==1
    %  [W,dW] = confArnoldi_inv(H,x);
    for k = 1:1

        [W,dW] = confArnoldi_inv(H,x);
        %      sk = -(W*d)./(dW*d);
        for j = 1:length(x)
            ff = W(j,:)*d;
            %            ff = f(x(j));
            df = dW(j,:)*d;

            sk = -ff./df;
            if abs(sk) > 1e-13
                x(j) = x(j) +   sk;
            end

        end
    end

end


end







% % % % % % % % % % % % % % % % % % % % % % % % % % % % % /
function [W,dW] = confArnoldi_inv(H,X)
% generate p(x) and p'(X) by  curculent coeffs
%
n = size(H,2);
X = X(:);
m = length(X);
W = zeros(m,n+1);          dW = W;
W(:,1) = X*0 + 1 ;

for k = 1:n
    w = X.*W(:,k);
    dw = W(:,k) + X.*dW(:,k);

    %          w = w - W(:,1:k)*H(1:k,k);
    %          W(:,k+1) = w/H(k+1,k);
    %
    %          dw = dw - dW(:,1:k)*H(1:k,k);
    %          dW(:,k+1) = dw/H(k+1,k);

    for j = 1:k
        w = w - W(:,j)*H(j,k);
        dw = dw - dW(:,j)*H(j,k);
    end
    W(:,k+1) = w/H(k+1,k);
    dW(:,k+1) = dw/H(k+1,k);
end

end