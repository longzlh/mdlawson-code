function [xi,ii] = zz_ratevalA(xout, a, b, H)
% Get the function handle of the rational approximant.
%   xi = p/q :        p = W*a; q = W*b;
% 
% cst =   norm(b),
% An = a(end)./b(end),

xout  = xout(:);

opt.clean = 1;


% -------------------
switch  opt.clean


    case 1
        q = zz_polyvalA(b,H, xout);
        p = zz_polyvalA(a,H, xout);
        xi = p./q;

        ii = find(abs(q)< 2e-14,1);

        if ~isempty(ii)
            warning('the output xi may not be  accurate, since q(xout) < eps! \n' );
        end


%  -----------------      Hopitol rule

% II = find(abs(q)< 1e-11);
% 
% II2 = find(abs(p(II)) < 1e-11) ; 
% 
% II = II(II2); 
% 
% Z = xout(II); 
% 
% [~,dW] = zz_confArnoldi_inv(H,Z); 
% 
% n1 = length(a);  n2 = length(b); 
% xi(II) = ( dW(:,1:n1)*a )./ ( dW(:,1:n2)*b); 

    
% -----------------  bad
    case 2

        n1 = find(abs(a)>6e-16, 1,'last');
        n2 = find(abs(b)>6e-16, 1,'last');
        H1 = H(1:n1,1:n1-1);
        H2 = H(1:n2,1:n2-1);


        % An = a(n1)./prod( diag(H1,-1) );
        % Bn = b(n2)./prod( diag(H2,-1) );

        if n1==n2
            AAA = a(n1)./b(n2); 
        elseif n1<n2
            AAA = a(n1)./b(n2)*prod(H2(n1+1:n2,n1:n2-1) );
        else    % n1>n2
            AAA = a(n1)./b(n2)./prod(H1(n2+1:n1,n2:n1-1) );
        end


        % n1 = length(a) - 1;   n2 = length(b) - 1;
        % q = zz_polyvalA(b,H(1:n2+1,1:n2), xout);
        % p = zz_polyvalA(a,H(1:n1+1,1:n1), xout);
        % % -------------------  delete common factor---------
        [zq, Aq] = zz_roots(H, b);
        [zp, Ap] = zz_roots(H, a);
        %
%         Ind = []; 
% for k = 1:length(zq)
%         zk = zq(k); 
% 
%         II = find(abs( xout-zk ) < 1e-11 ); 
%         Ind = [Ind; II]; 
% end

 
        mj = zp - zq.';



        [J,K] = find(abs(mj) < 1e-12);

        Xp = zp;
        Xq = zq;
        %
        if ~isempty(J)
            fprintf('the common factors are deleted ! \n' );
%             Xp(J) - Xq(K)
            Xp(J) = [];
            Xq(K) = [];
        end
        p = prod(xout(:) - Xp.', 2);
        q = prod(xout(:) - Xq.', 2);


        xi = AAA*(p./q);      ii = [];
        %  ------------------------------------------------------------------



end

