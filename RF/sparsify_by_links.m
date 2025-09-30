function K_sparse = sparsify_by_links(K_dense, num_links)

    K_sparse = zeros(size(K_dense));
    [~, idx] = sort(abs(K_dense(:)), 'descend');
    keep = idx(1:min(num_links, numel(idx)));
    K_sparse(keep) = K_dense(keep);
end
