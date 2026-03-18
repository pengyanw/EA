function L = laplacian_matrix(g)
%LAPLACIAN_MATRIX 计算 graph 类型 g 的拉普拉斯矩阵，不使用内置 laplacian()

    A = adjacency(g, 'weighted');   % 邻接矩阵
    D = diag(sum(A, 2));            % 度矩阵
    L = D - A;                      % 拉普拉斯矩阵
end
