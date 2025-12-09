function [K_best, history] = post_refine_K(A, B, Q, R, Adj, K_init, maxIter)
%--------------------------------------------------------------------------
% Post-Processing refinement of EA's final controller.
% Inputs:
%   A, B         system matrices
%   Q, R         weight matrices of the individual
%   Adj          adjacency matrix
%   K_init       final sparse K from EA
%   maxIter      max refinement iterations
%
% Outputs:
%   K_best       refined K
%   history      cost history
%--------------------------------------------------------------------------

if nargin < 7
    maxIter = 100;
end

Nx = size(A,1);
Nu = size(B,2);

% shortest-path distances
G = graph(Adj);
D = distances(G);

% reference dense LQR
K_dense = -dlqr(A,B,Q,R);

% cost using Frobenius distance
cost_fun = @(K) norm(K - K_dense, 'fro');

K_best = K_init;
best_cost = cost_fun(K_best);
history = zeros(maxIter,1);

num_links = nnz(K_init);

for t = 1:maxIter
    
    % perform one graph-driven local swap
    K_new = graph_local_swap(A, B, Q, R, Adj, K_best, num_links, D);
    
    % check stability
    if max(abs(eig(A + B*K_new))) >= 1
        history(t) = best_cost;
        continue;
    end
    
    % evaluate new cost
    new_cost = cost_fun(K_new);
    history(t) = best_cost;
    
    % accept if improved
    if new_cost < best_cost
        best_cost = new_cost;
        K_best = K_new;
    end
end

end
