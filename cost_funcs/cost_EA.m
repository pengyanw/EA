function J_total = cost_EA(A, B, Q, R, K_sparse, costBM, alpha, ~, gers_weight)
% COST_EA  Compute the evolutionary algorithm cost for a given controller.
%
% Usage:
%   J_total = cost_EA(A, B, Q, R, K_sparse)
%   J_total = cost_EA(A, B, Q, R, K_sparse, costBM, alpha, max_links)
%   J_total = cost_EA(A, B, Q, R, K_sparse, costBM, alpha, max_links, gers_weight)
%
% Inputs:
%   A, B        - System matrices
%   Q, R        - LQR weight matrices (default diagonal)
%   K_sparse    - Controller gain matrix (Nu x Nx), may contain zeros
%
% Optional:
%   costBM      - Benchmark dense LQR cost (default = 1)
%   alpha       - Weighting factor between LQR and hardware cost (default = 0.5)
%   max_links   - Normalization term for nnz(K) (default = numel(K_sparse))
%   gers_weight - Gershgorin soft penalty weight (default = 0, disabled)
%                 Typical range: 1~10. Penalizes Gershgorin radius exceeding 1.
%
% Output:
%   J_total     - Normalized total cost value
%
% Components of cost:
%   1) Normalized LQR performance cost
%   2) Communication penalty (nnz(K))
%   3) Actuator penalty (nonzero rows)
%   4) Sensor penalty (nonzero cols)
%   5) Gershgorin stability margin penalty (if gers_weight > 0)
%
% Example:
%   J = cost_EA(A,B,Q,R,K_sparse, costBM, 0.6, Nu*Nx);
%   J = cost_EA(A,B,Q,R,K_sparse, costBM, 0.6, Nu*Nx, 5);  % with Gershgorin
%
% Author: Pengyang Wu (EA-LQR framework)
% ================================================================

if nargin < 9, gers_weight = 0; end

% ----------- Step 1: Closed-loop stability check -----------
eigVals = eig(A + B*K_sparse);
if max(abs(eigVals)) >= 1.0
    J_total = 1e9;
    return;
end

% ----------- Step 2: Compute normalized LQR performance cost -----------
try
    lqr_perf_cost = get_lqr_cost(A, B, Q, R, K_sparse);
catch
    J_total = 1e12;
    return;
end
J_perf_norm = lqr_perf_cost / costBM;

% ----------- Step 3: Structural sparsity penalties -----------
comm_cost   = nnz(K_sparse);
active_rows = nnz(any(K_sparse,2));
active_cols = nnz(any(K_sparse,1));

% ----------- Step 4: Weighted total cost -----------
w_comm = 0.05*(1 - alpha);
w_row  = 0.4*(1 - alpha);
w_col  = 0.2*(1 - alpha);

J_total =  J_perf_norm + ...
           w_comm * comm_cost + ...
           w_row  * active_rows + ...
           w_col  * active_cols;

% ----------- Step 5: Gershgorin stability margin penalty -----------
% if gers_weight > 0
%     L_gers  = gersgorin_penalty(A, B, K_sparse);
%     J_total = J_total + gers_weight * L_gers;
% end

end
