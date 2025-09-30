function models = rf_train_K_models(gene_mat, K_cell, Nu, Nx)
% Train RF regressors for each K(i,j).
% gene_mat: [Nsample x geneLength], every row is a combination of [diag(Q), diag(R), num_links]
% K_cell  : {Nsample x 1}, every element is: Nu x Nx K matrix
% Nu, Nx   
%
% models  : Nu x Nx  cell, models{i,j} RF

    N = size(gene_mat,1);
    assert(numel(K_cell)==N, 'gene_mat, K_cell mismatch');

    
    Y = cell(Nu, Nx);
    for i = 1:Nu
        for j = 1:Nx
            Y{i,j} = zeros(N,1);
        end
    end
    for s = 1:N
        K = K_cell{s};
        for i = 1:Nu
            for j = 1:Nx
                Y{i,j}(s) = K(i,j);
            end
        end
    end

    % train every (i,j): RF{i,j}
    t = templateTree('MinLeafSize',5,'MaxNumSplits',min(50, size(gene_mat,1)-1));
    models = cell(Nu, Nx);
    for i = 1:Nu
        for j = 1:Nx
            models{i,j} = fitrensemble(gene_mat, Y{i,j}, ...
                'Method','Bag', 'NumLearningCycles',100, 'Learners',t);
        end
    end
end
