function D = util_matcell(M,C)
%MATCELL    Apply linear transformation matrix M to a column cell 
%           array of coefficient matrices C.

O = 0*C{1}; % sparse or dense
D = cell(size(M,1),1);
for i = 1:size(M,1),
    D{i} = O;
    for j = 1:size(M,2),
        D{i} = D{i} + M(i,j)*C{j};
    end
end