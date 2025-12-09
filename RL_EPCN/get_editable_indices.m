function idx = get_editable_indices(K_dense, adj)
    % 例如选择 |K_dense| 最大的前 50 个
    [~,ord] = sort(abs(K_dense(:)),'descend');
    top = ord(1:50);
    [I,J]=ind2sub(size(K_dense),top);
    idx = [I,J];
end
