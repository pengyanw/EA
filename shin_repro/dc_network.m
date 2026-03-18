function G = dc_network()
    mpc = pglib_opf_case118_ieee();  % ensure this function is on path

    c = (132e3)^2 / 1e5 * (1e-9);  % scale

    n = size(mpc.bus, 1);

    % Use containers.Map to keep only the last-seen susceptance per node pair
    % (matches Julia SimpleWeightedGraph which overwrites duplicate edges)
    edge_sus = containers.Map('KeyType','char','ValueType','double');
    for i = 1:size(mpc.branch, 1)
        from = mpc.branch(i, 1);
        to   = mpc.branch(i, 2);
        % column 3 = BR_R (resistance), column 4 = BR_X (reactance)
        % susceptance = Im(1/(R+jX))
        R_ij = mpc.branch(i, 3);
        X_ij = mpc.branch(i, 4);
        b = abs(imag(1 / (R_ij + 1j * X_ij)));  % branch susceptance
        key = sprintf('%d_%d', min(from,to), max(from,to));
        edge_sus(key) = b * c;   % overwrites duplicate — matches Julia
    end

    % Build simple graph from deduplicated edges
    G = graph();
    G = addnode(G, n);
    ks = keys(edge_sus);
    for ki = 1:length(ks)
        parts = sscanf(ks{ki}, '%d_%d');
        G = addedge(G, parts(1), parts(2), edge_sus(ks{ki}));
    end
end
