function blocksizes = util_pencil_blocksizes(K,H,s0)
%UTIL_PENCIL_BLOCKSIZES    Infer blocksizes in a Hessenberg matrix pencil in thin form.
%
% Identify the block sizes of a block upper-Hessenberg pencil (H, K), 
% where s0 is the dimension of the first diagonal block, i.e., the number 
% of starting vectors.
%
% Assumes Hessenberg matrix pencil is in fat form.

blocksizes = [1,s0+1];
i = 2;
[n,m] = size(K);
col = [1:s0:m+1];
while blocksizes(end) < size(K,1)
    
    if 2*blocksizes(i) <= n
        l = 2*blocksizes(i)-blocksizes(i-1)-1;
    else
        l = n;
    end
    
    if K(blocksizes(i):l, col(i-1):col(i)-1) == zeros(l-blocksizes(i)+1, col(i)-col(i-1))% if pole = inf
        c_u = blocksizes(i) + rank(H(blocksizes(i):l, col(i-1):col(i)-1));
    else
        c_u = blocksizes(i) + rank(K(blocksizes(i):l, col(i-1):col(i)-1)); % All other cases
    end
    
    blocksizes = [blocksizes,c_u];
    i = i+1;
end

if blocksizes(i) <= size(K,1)
    blocksizes(i+1) = size(K,1)+1;
end

blocksizes = diff(blocksizes);
end