function [R,rmse,varargout] = fixed_support_block_aaa(F,pts,Z,opts)
%FIXED_SUPPORT_BLOCK_AAA Block-AAA barycentric fit with prescribed supports.
%
% [R,rmse,poles,zeros,out] = fixed_support_block_aaa(F,pts,Z,opts) computes the
% block-AAA barycentric rational approximant using the prescribed support
% points Z. The points in Z must be selected from pts. This is the fixed
% support-point part of block_aaa, with the adaptive point-selection loop
% removed.
%
% For backward compatibility, [R,rmse,out] returns the output structure as
% the third output. Poles and zeros are always stored in out.poles and
% out.zeros.
%
% Inputs:
%   F    function handle @(z) returning an m-by-n matrix, or a cell array
%        of length(pts) containing the sampled m-by-n matrices.
%   pts  distinct sampling points.
%   Z    prescribed support points, selected from pts.
%   opts optional structure:
%        opts.matchTol -- tolerance for matching Z to pts (default 100*eps).
%
% Outputs:
%   R    function handle with
%        R(z) = (sum_k Dk/(z-Zk)) \ (sum_k Ck/(z-Zk)),
%        Ck = Dk*F(Zk).
%   rmse root mean squared error on all pts.
%   poles poles of the barycentric denominator.
%   zeros zeros/nonlinear eigenvalues of R. For non-square R, this is [].
%   out  structure containing zk, supportIndex, Ck, Dk, errors, poles, zeros.

if nargin < 4
    opts = struct();
end
if ~isfield(opts,'matchTol')
    opts.matchTol = 100*eps;
end

pts = pts(:);
Z = Z(:);
npts = length(pts);
numSupports = length(Z);

if npts == 0
    error('fixed_support_block_aaa:emptyPts', 'pts must not be empty.');
end
if numSupports == 0
    error('fixed_support_block_aaa:emptySupports', 'Z must not be empty.');
end
if numSupports >= npts
    error('fixed_support_block_aaa:tooManySupports', ...
        'Z must contain fewer points than pts so that least-squares points remain.');
end

if iscell(F)
    if numel(F) ~= npts
        error('fixed_support_block_aaa:badCellData', ...
            'When F is a cell array, numel(F) must equal length(pts).');
    end
    FF = F(:);
else
    FF = cell(npts,1);
    for i = 1:npts
        FF{i} = F(pts(i));
    end
end

[m,n] = size(FF{1});
for i = 2:npts
    if ~isequal(size(FF{i}), [m,n])
        error('fixed_support_block_aaa:inconsistentSize', ...
            'All sampled matrices F(pts(i)) must have the same size.');
    end
end

zk_ind = local_match_support_indices(pts,Z,opts.matchTol);
lamind = true(npts,1);
lamind(zk_ind) = false;
lamind = find(lamind).';

zk = pts(zk_ind);
lam = pts(lamind);
p = length(lam);

M = zeros(numSupports*m,p*n);
for j = 1:numSupports
    rowRange = (j-1)*m+1:j*m;
    for i = 1:p
        colRange = (i-1)*n+1:i*n;
        M(rowRange,colRange) = (FF{lamind(i)} - FF{zk_ind(j)}) ...
            /(lam(i) - zk(j));
    end
end

EF = local_left_null_basis(M,m);

Ck = cell(numSupports,1);
Dk = cell(numSupports,1);
for j = 1:numSupports
    Dk{j} = EF(:,(j-1)*m+1:j*m);
    Ck{j} = Dk{j}*FF{zk_ind(j)};
end

R = @(z) local_eval_bary(z,zk,Ck,Dk);
poleValues = local_nonlinear_eig(Dk,zk);
if m == n
    zeroValues = local_nonlinear_eig(Ck,zk);
    zeroStatus = 'computed';
else
    zeroValues = [];
    zeroStatus = 'not computed: zeros are defined here only for square matrix functions';
end

pointErrors = zeros(npts,1);
for i = 1:npts
    pointErrors(i) = norm(R(pts(i)) - FF{i},'fro')^2;
end
pointErrors(~isfinite(pointErrors)) = 0;
rmse = sqrt(sum(pointErrors)/npts);

out = struct();
out.zk = zk;
out.supportIndex = zk_ind;
out.lsIndex = lamind;
out.Ck = Ck;
out.Dk = Dk;
out.loewnerMatrix = M;
out.pointErrors = pointErrors;
out.supportRmse = sqrt(sum(pointErrors(zk_ind))/numSupports);
out.poles = poleValues;
out.zeros = zeroValues;
out.zeroStatus = zeroStatus;

if nargout == 3
    varargout{1} = out;
elseif nargout >= 4
    varargout{1} = poleValues;
    varargout{2} = zeroValues;
    if nargout >= 5
        varargout{3} = out;
    end
end

end

function supportIndex = local_match_support_indices(pts,Z,matchTol)
supportIndex = zeros(length(Z),1);
for j = 1:length(Z)
    scale = max(1,abs(Z(j)));
    [dist,idx] = min(abs(pts - Z(j)));
    if dist > matchTol*scale
        error('fixed_support_block_aaa:supportNotInPts', ...
            'Every support point in Z must be selected from pts.');
    end
    supportIndex(j) = idx;
end
if numel(unique(supportIndex)) ~= numel(supportIndex)
    error('fixed_support_block_aaa:duplicateSupports', ...
        'Z must contain distinct support points.');
end
end

function EF = local_left_null_basis(M,m)
if size(M,2) >= size(M,1)
    [~,~,V] = svd(M',0);
else
    [~,~,V] = svd(M');
end
if size(V,2) < m
    error('fixed_support_block_aaa:notEnoughBasisVectors', ...
        'Could not compute enough barycentric denominator weight vectors.');
end
EF = V(:,end-m+1:end)';
end

function Rz = local_eval_bary(z,zk,Ck,Dk)
N = zeros(size(Ck{1}));
D = zeros(size(Dk{1}));

[val,ind] = min(abs(z-zk));
if val < 10*eps
    Rz = Dk{ind}\Ck{ind};
    return
end

for j = 1:length(zk)
    N = N + Ck{j}/(z-zk(j));
    D = D + Dk{j}/(z-zk(j));
end
Rz = D\N;
end

function evs = local_nonlinear_eig(C,z)
[m,n] = size(C{1});
if m ~= n
    error('fixed_support_block_aaa:nonSquareEigenProblem', ...
        'Pole/zero finding works only for square matrix coefficients.');
end

ell = length(z)-1;
if ell <= 0
    evs = [];
    return
end

Cscaled = C;
w = zeros(ell+1,1);
for i = 1:ell+1
    x = z - z(i);
    x(i) = [];
    w(i) = 1/prod(x);
    Cscaled{i} = Cscaled{i}/w(i);
end

theta = w(1:end-1)./w(2:end);

A0 = zeros(ell-1,ell);
for i = 1:ell-1
    A0(i,i) = z(i);
    A0(i,i+1) = -z(i+2)*theta(i);
end
L0bott = kron(A0,eye(m,m));

L0top = [];
for i = 1:ell-1
   L0top = [L0top, z(i+1)*Cscaled{i}]; %#ok<AGROW>
end
L0top = [L0top, z(ell+1)*Cscaled{ell} ...
    + z(ell)*(1/theta(ell))*Cscaled{ell+1}];
L0 = [L0top; L0bott];

A1 = zeros(ell-1,ell);
for i = 1:ell-1
    A1(i,i) = 1;
    A1(i,i+1) = -theta(i);
end
L1bott = kron(A1,eye(m,m));

L1top = [];
for i = 1:ell-1
   L1top = [L1top, Cscaled{i}]; %#ok<AGROW>
end
L1top = [L1top, Cscaled{ell} + (1/theta(ell))*Cscaled{ell+1}];
L1 = [L1top; L1bott];

evs = sort(eig(L0,L1));
end
