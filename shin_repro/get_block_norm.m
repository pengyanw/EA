function block_norm = get_block_norm(K, G, uloc)
% uloc = control indices (K rows)
% block_norm[k+1] stores max ||K(i,j)|| at distance k

    D = distances(G);  % all-pairs distances
    max_dist = max(D(:));
    block_norm = zeros(1, max_dist + 1);

    for ii = 1:length(uloc)
        i = uloc(ii);          % control index
        d = D(i, :);           % distances from control i to all states
        for j = 1:size(K, 2)   % loop over states
            dist = d(j);
            block_norm(dist + 1) = max(block_norm(dist + 1), abs(K(i,j)));
        end
    end
end
