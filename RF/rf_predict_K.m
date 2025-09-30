function K_pred = rf_predict_K(models, gene_row)

% models : Nu x Nx  cell (rf_train_K_models )
% gene_row: [1 x geneLength]
% K_pred : Nu x Nx

    [Nu, Nx] = size(models);
    K_pred = zeros(Nu, Nx);
    for i = 1:Nu
        for j = 1:Nx
            K_pred(i,j) = predict(models{i,j}, gene_row);
        end
    end
end
