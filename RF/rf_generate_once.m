function [K_new, J_new, isStable, specRad] = rf_generate_once( ...
    models, gene_row, A, B, Q_bm, R_bm, costBM, alpha, max_links)

    
    % assume gene = [diag(Q), diag(R), num_links]
    num_links = round(gene_row(end));

    % predict dense K and sparsify
    K_dense  = rf_predict_K(models, gene_row);
    K_sparse = sparsify_by_links(K_dense, num_links);

    % loss
    [J_new, isStable, specRad] = eval_K_cost(K_sparse, A, B, Q_bm, R_bm, costBM, alpha, max_links);
    K_new = K_sparse;
end
