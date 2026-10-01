function check = test_fixed_support_block_aaa()
%TEST_FIXED_SUPPORT_BLOCK_AAA Basic checks for prescribed-support block-AAA.

pts = 1i*linspace(0.1,1,40).';
F = @(z) [1/(z+1),      1/(z^2-5),                 z/(z+3); ...
          1/(z+2+1i),  (2+z^2)/(z^3+3*z^2+1),     7/(z+4)];

supportIndex = [1, 2, 8, 16, 27, 40];
Z = pts(supportIndex);

[R,fitRmse,poleValues,zeroValues,out] = fixed_support_block_aaa(F,pts,Z);

FF = cell(length(pts),1);
for i = 1:length(pts)
    FF{i} = F(pts(i));
end
[Rcell,cellRmse,outCell] = fixed_support_block_aaa(FF,pts,Z);

interpErr = 0;
for j = 1:length(Z)
    interpErr = max(interpErr,norm(R(Z(j)) - F(Z(j)),'fro'));
end

recomputedRmse = local_rmse(pts,R,FF);
cellDiff = 0;
for i = 1:length(pts)
    cellDiff = max(cellDiff,norm(R(pts(i)) - Rcell(pts(i)),'fro'));
end

check = interpErr < 1e-10 ...
    && abs(fitRmse - recomputedRmse) < 1e-14 ...
    && abs(fitRmse - cellRmse) < 1e-14 ...
    && cellDiff < 1e-12 ...
    && isequal(poleValues,out.poles) ...
    && isempty(zeroValues) ...
    && isempty(out.zeros) ...
    && isequal(out.supportIndex(:),supportIndex(:)) ...
    && isequal(outCell.supportIndex(:),supportIndex(:));

Fsquare = @(z) [1/(z+1), z/(z+2); 1/(z+3), (z+1)/(z+4)];
[~,~,squarePoles,squareZeros,squareOut] = fixed_support_block_aaa(Fsquare,pts,Z);
check = check ...
    && ~isempty(squarePoles) ...
    && ~isempty(squareZeros) ...
    && isequal(squarePoles,squareOut.poles) ...
    && isequal(squareZeros,squareOut.zeros);

fprintf('fixed_support_block_aaa: rmse = %.4e, support interpolation err = %.4e\n', ...
    fitRmse, interpErr);

end

function err = local_rmse(pts,F1,F2)
npts = length(pts);
FF1 = local_sample_values(pts,F1);
FF2 = local_sample_values(pts,F2);

err = 0;
for i = 1:npts
    err = err + norm(FF1{i} - FF2{i},'fro')^2;
end
err = sqrt(err/npts);
end

function FF = local_sample_values(pts,F)
npts = length(pts);
if iscell(F)
    FF = F(:);
else
    FF = cell(npts,1);
    for i = 1:npts
        FF{i} = F(pts(i));
    end
end
end
