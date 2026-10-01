function R = feval(obj, varargin)
%FEVAL    Evaluate RKFUNB at scalar or or matrix arguments.
%
% Calling syntax: 
%   - R = feval(obj, z)  -- evaluate at scalar z
%   - R = feval(obj, Z)  -- evaluate for an sxs matrix Z (s = block size)
%   - R = feval(obj,Z,B) -- evaluate for a matrix times vector (circ op)

s = size(obj.coeffs{1},1); % block size

if nargin == 2,
    A = varargin{1};
    if ~isnumeric(A),
        error('FEVAL: RKFUNB evaluation only possible for numeric argument.');
    end
    if norm(size(A) - [1,1],inf) == 0, % scalar A?
        A = A*eye(s);
    end
    if norm(size(A) - [s,s],inf) > 0, % check A square
        error('FEVAL: RKFUNB evaluation only possible when single argument is square and of correct size.');
    end
    B = eye(s);
end

if nargin == 3,
    A = varargin{1};
    B = varargin{2};
    if ~isnumeric(A) || ~isnumeric(B),
        error('FEVAL: RKFUNB evaluation only possible for numeric arguments.');
    end
    if norm(size(A) - [1,1],inf) == 0, % scalar A
        A = A*speye(size(B,1));
    end
    if diff(size(A)) > 0 || size(A,1) ~= size(B,1), % check A square and compatible with B
        error('FEVAL: RKFUNB evaluation only possible when first argument is square and compatible with second argument.');
    end
    if size(B,2) ~= s, % check B is of correct width
        error('FEVAL: RKFUNB evaluation only possible when second argument is of correct block size.');
    end
end



[V,K,H] = rat_krylov(A,B,obj.K,obj.H);
%norm(A*V*K - V*H)

R = 0;
for j = 1:round(size(V,2)/s), % m+1
   R = R + V(:,1+(j-1)*s:j*s)*obj.coeffs{j};
end


end

