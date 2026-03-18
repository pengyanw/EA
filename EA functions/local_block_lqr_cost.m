function [Jloc, blocks] = local_block_lqr_cost(adjMtx, A, B, Q, R, K, W)
% LOCAL_BLOCK_LQR_COST
% Blockwise approximation of LQR cost using graph-distance-based neighborhoods
%
% INPUTS:
%   adjMtx : n x n adjacency matrix
%   A,B,Q,R,K : system matrices and controller
%   W      : disturbance covariance (optional, default = I)
%
% OUTPUTS:
%   Jloc   : approximated LQR cost
%   blocks : cell array of local state index sets S_i

    if nargin < 7 || isempty(W)
        W = eye(size(A,1));
    end

    n = size(A,1);

    %% 1. Graph distance
    G = graph(adjMtx);
    D = distances(G);        % n x n shortest-path distances
    D(isinf(D)) = inf;

    %% 2. Closed-loop matrices
    Acl = A + B*K;
    Qc  = Q + K'*R*K;

    %% 3. Build local blocks
    blocks = cell(n,1);
    weights = cell(n,1);

    for i = 1:n
        % graph distances from node i
        di = D(i,:);

        % decay weights
        wi = exp(-di);

        % ignore self
        wi(i) = 0;

        % number of neighbors to keep: half of degree
        deg_i = nnz(adjMtx(i,:));
        keepN = max(1, floor(deg_i/2));

        % pick strongest connections (smallest distance)
        [~, idx] = sort(wi, 'descend');
        Si = idx(1:keepN);

        % include center node
        Si = unique([i, Si]);

        blocks{i}  = Si;
        weights{i} = wi(Si);
    end

    %% 4. Coverage count for normalization
    cover_count = zeros(n,1);
    for i = 1:n
        cover_count(blocks{i}) = cover_count(blocks{i}) + 1;
    end
    inv_cover = 1 ./ max(cover_count,1);

    %% 5. Local Lyapunov solves (can be parfor)
    contrib = cell(n,1);

    parfor i = 1:n
        Si = blocks{i};

        Ai = Acl(Si,Si);
        Qi = Qc(Si,Si);
        Wi = W(Si,Si);

        % optional damping to avoid optimistic bias
        % eta = 0.1;
        % Ai = Ai - eta*eye(length(Si));

        Pi = lyap(Ai', Qi);

        % local diagonal contribution
        contrib{i} = diag(Pi * Wi);
    end

    %% 6. Aggregate (method 2)
    J_state = zeros(n,1);
    for i = 1:n
        Si = blocks{i};
        J_state(Si) = J_state(Si) + contrib{i};
    end

    Jloc = sum(J_state .* inv_cover);

end
