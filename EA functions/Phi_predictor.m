function [Jhat, out] = Phi_predictor(A, B, K_dense, adjMtx, linkTraj, J0, w_c, ea_params, opts)
% PHI_PREDICTOR  Per-generation EA cost prediction from the drift bound of Theorem 1 (eq:plotted_bound).
%
%   This implements the convergence result of Section IV directly, rather than
%   fitting an exponential envelope to the EA curve it is supposed to predict.
%   Every constant below is estimated from (A, B, K_dense) and the graph, so the
%   resulting curve is falsifiable.
%
%   Definitions implemented (equation numbers refer to the final paper):
%
%     Lemma 1     |K_d,ij| <= Upsilon * rho^{d_G(i,j)}
%     (10)        Phi(h)   = (L_J * Upsilon^2 * rho^{2h} / J_LQR(K_d))
%                            * (sqrt(Nu*Nx) + 1/2)
%     (11)        h*       = (1/(2|log rho|))
%                            * log( L_J*Upsilon^2*(sqrt(Nu*Nx)+1/2)
%                                   / (w_c*J_LQR(K_d)) )
%     Def. 2(c)   h_t      = largest integer such that every node pair at graph
%                            distance <= h_t has an edge in the controller
%                            adjacency matrix of K_s([ell_t, 1, 1])
%     (17)        p_imp(t) = 1{h_{t-1} > h*} / (Np * (2d+1))
%                 P_imp(t) = 1 - (1 - p_imp(t))^(Np - ne)
%     Thm. 1     E[J_EA(t-1) - J_EA(t)] >= (w_c - Phi(h_{t-1})) * P_imp(t)
%
%   The recursion is seeded at the measured J_EA of generation 1 and driven by
%   the MEASURED link-count trajectory: the bound lower-bounds the per-generation
%   drift given the current state, it does not predict ell_t itself.
%
% Syntax:
%   [Jhat, out] = Phi_predictor(A, B, K_dense, adjMtx, linkTraj, J0, w_c, ea_params)
%   [Jhat, out] = Phi_predictor(..., opts)
%
% Inputs:
%   A, B        - Plant matrices (Nx x Nx, Nx x Nu)
%   K_dense     - Dense LQR gain K_d (Nu x Nx)
%   adjMtx      - Node adjacency matrix of the system graph (numNodes x numNodes)
%   linkTraj    - Per-generation best link count ell*_t (maxGen x 1), from
%                 history.bestLinks
%   J0          - Measured J_EA of the best individual at generation 1
%                 (raw, unnormalized units)
%   w_c         - Communication weight (0.05 in the paper)
%   ea_params   - Struct with popSize (Np), nTop (ne), and optionally mutStep (d)
%   opts        - (optional) struct:
%                 .Q, .R       - LQR weights (default identity)
%                 .d           - Mutation step bound (default 5)
%                 .statesPerNode - States per node (default 2: theta_i, omega_i)
%                 .verbose     - Print the estimated constants (default true)
%
% Outputs:
%   Jhat  - Predicted cost trajectory (maxGen x 1), raw units. Because the bound is a
%           lower bound on the per-generation DECREASE, Jhat is an upper envelope
%           of the true best-cost curve: Jhat(t) >= J_EA(theta*_t) for all t.
%   out   - Struct with fields Upsilon, rho, L_J, h_star, h_traj, Phi_traj,
%           P_imp, J_LQR_Kd, looseness.
%
% Author: Pengyang Wu

if nargin < 9, opts = struct(); end
Q             = get_opt(opts, 'Q',             eye(size(A,1)));
R             = get_opt(opts, 'R',             eye(size(B,2)));
d             = get_opt(opts, 'd',             5);
statesPerNode = get_opt(opts, 'statesPerNode', 2);
verbose       = get_opt(opts, 'verbose',       true);

[Nu, Nx]  = size(K_dense);
numNodes  = size(adjMtx, 1);
maxGen    = numel(linkTraj);

Np = ea_params.popSize;
ne = ea_params.nTop;

%% ---------- Graph distances between nodes ----------
G = graph(adjMtx ~= 0);
D = distances(G);                       % numNodes x numNodes hop counts

% Controller entry K(u, x): row u is the actuator at node u, column x belongs to
% node colNode(x). The two plant families in this repo use DIFFERENT state
% orderings, so this must not be hard-coded:
%   grid      : interleaved, node ii owns states [2*ii-1, 2*ii]  -> ceil(x/2)
%   IEEE 13-bus: blocked,    node ii owns states [ii, N+ii]      -> mod-based
% Pass opts.colNode explicitly for anything that is not the interleaved layout.
rowNode = (1:Nu)';
colNode = get_opt(opts, 'colNode', ceil((1:Nx) / statesPerNode));
if numel(colNode) ~= Nx
    error('Phi_predictor:colNode', ...
        'opts.colNode must have %d entries (one per state), got %d.', Nx, numel(colNode));
end
Dk      = D(rowNode, colNode);          % Nu x Nx entrywise hop distance

%% ---------- (a) Estimate Upsilon and rho from Lemma 1 ----------
% Lemma 1 is an upper envelope, so regress the per-distance MAXIMUM magnitude:
%   log max_{d_G(i,j)=h} |K_d,ij|  =  log Upsilon + h * log rho
hVals = unique(Dk(isfinite(Dk)))';
mMax  = nan(size(hVals));
for k = 1:numel(hVals)
    v = abs(K_dense(Dk == hVals(k)));
    v = v(v > 0);
    if ~isempty(v), mMax(k) = max(v); end
end
ok = isfinite(mMax) & mMax > 0;
if nnz(ok) < 2
    error('Phi_predictor:decayFit', 'Not enough distinct hop distances to fit Lemma 1.');
end
coef    = polyfit(hVals(ok), log(mMax(ok)), 1);
rho     = exp(coef(1));
Upsilon = exp(coef(2));

% Lemma 1 requires rho in (0,1) and Upsilon >= 1.  Clamp and warn rather than
% silently producing a bound that the lemma does not support.
if rho >= 1 || rho <= 0
    warning('Phi_predictor:rho', ...
        'Fitted rho = %.4f is outside (0,1); Lemma 1 does not hold on this fit.', rho);
    rho = min(max(rho, 1e-6), 1 - 1e-6);
end
if Upsilon < 1
    % Rescale to the smallest envelope with Upsilon >= 1 that still dominates the data.
    Upsilon = 1;
end
% Make the envelope valid by construction: inflate Upsilon until it dominates
% every entry, otherwise Phi(h) is not an upper bound on the pruning cost.
slack   = max(abs(K_dense(:)) ./ (Upsilon * rho.^Dk(:)));
if slack > 1, Upsilon = Upsilon * slack; end

%% ---------- (b) Estimate L_J along the pruning path ----------
% Eq. (8) is  J_LQR(K) - J_LQR(K_d) <= (L_J/2)*||K - K_d||_F^2, so the smallest
% constant that makes (8) hold on the family the EA actually visits is
%   L_J = max_ell  2*(J_LQR(Pi_ell(K_d)) - J_LQR(K_d)) / ||Pi_ell(K_d) - K_d||_F^2.
%
% Lemma 2 states the Lipschitz-gradient property only on a SUBLEVEL SET
% S = {K : J_LQR(K) <= c}, and that restriction is essential: J_LQR blows up at
% the stability boundary, so sampling the whole pruning path yields a
% meaningless L_J (~1e30).  The relevant sublevel set here is the one the elite
% actually occupies: every term of J_EA is non-negative, so along the elite
% sequence  J_LQR(K)/J_LQR(K_d) <= J_EA(theta*_0) = J0.
J_LQR_Kd = get_lqr_cost(A, B, Q, R, K_dense);
c_sub    = get_opt(opts, 'sublevel', J0) * J_LQR_Kd;

[~, sortIdx] = sort(abs(K_dense(:)), 'descend');
nnzKd   = nnz(K_dense);
ellGrid = unique(round(linspace(1, nnzKd, min(120, nnzKd))));

L_J     = 0;
nInSub  = 0;
for ell = ellGrid
    K_s = zeros(Nu, Nx);
    keep = sortIdx(1:ell);
    K_s(keep) = K_dense(keep);
    dK2 = norm(K_s - K_dense, 'fro')^2;
    if dK2 < eps, continue; end
    Jl = get_lqr_cost(A, B, Q, R, K_s);
    if ~isfinite(Jl) || Jl > c_sub, continue; end   % outside S: (8) says nothing
    nInSub = nInSub + 1;
    L_J = max(L_J, 2 * (Jl - J_LQR_Kd) / dK2);
end
if L_J <= 0
    warning('Phi_predictor:LJ', ...
        'No sampled pruned controller lay in the sublevel set; L_J fell back to 1.');
    L_J = 1;
end

%% ---------- (c) Effective truncation distance h_t from ell*_t ----------
% Definition 2(c): build K_s([ell, 1, 1]) = top-ell entries of K_d by magnitude,
% form its controller adjacency matrix over nodes, then h is the largest integer
% such that EVERY node pair within h hops carries an edge.
maxD    = max(D(isfinite(D)));
h_of_ell = containers.Map('KeyType', 'double', 'ValueType', 'double');
h_traj  = zeros(maxGen, 1);
for t = 1:maxGen
    ell = max(1, min(round(linkTraj(t)), nnzKd));
    if isKey(h_of_ell, ell)
        h_traj(t) = h_of_ell(ell);
        continue;
    end
    K_s = false(Nu, Nx);
    K_s(sortIdx(1:ell)) = true;

    % Node-level controller adjacency: edge (i,j) iff block K[row i, cols j] ~= 0.
    % Only the Nu actuator rows can carry edges, so restrict the comparison to
    % those rows; pairs (i,j) with i beyond Nu are not part of E_K by definition.
    E  = accumarray_blocks(K_s, colNode, numNodes) > 0;   % Nu x numNodes
    Dn = D(rowNode, :);                                   % Nu x numNodes

    missing = Dn(~E & isfinite(Dn));
    if isempty(missing)
        hh = maxD;
    else
        hh = min(missing) - 1;   % -1 means even a 0-hop (self) block is absent
    end
    h_of_ell(ell) = hh;
    h_traj(t)     = hh;
end

%% ---------- (d) Phi, h*, and the Theorem 1 recursion ----------
scale    = sqrt(Nu * Nx) + 0.5;
Phi      = @(h) (L_J * Upsilon^2 * rho.^(2*h) / J_LQR_Kd) * scale;
h_star   = log(L_J * Upsilon^2 * scale / (w_c * J_LQR_Kd)) / (2 * abs(log(rho)));

Phi_traj = Phi(h_traj);
p_imp    = (1 / (Np * (2*d + 1))) * double(h_traj > h_star);
P_imp    = 1 - (1 - p_imp).^(Np - ne);

Jhat    = zeros(maxGen, 1);
Jhat(1) = J0;
for t = 2:maxGen
    drift   = max(w_c - Phi_traj(t-1), 0) * P_imp(t-1);
    Jhat(t) = Jhat(t-1) - drift;
end

% How small would L_J have to be for the certified region {h_t > h*} to be
% non-empty at all?  h* < max(h_traj) requires
%   L_J < w_c*J_LQR(K_d)*rho^{-2*max(h_traj)} / (Upsilon^2 * scale).
L_J_crit = w_c * J_LQR_Kd * rho^(2*max(h_traj)) / (Upsilon^2 * scale);

out = struct('Upsilon', Upsilon, 'rho', rho, 'L_J', L_J, 'h_star', h_star, ...
             'h_traj', h_traj, 'Phi_traj', Phi_traj, 'p_imp', p_imp, ...
             'P_imp', P_imp, 'J_LQR_Kd', J_LQR_Kd, 'L_J_crit', L_J_crit, ...
             'nCertifiedGens', nnz(h_traj > h_star), 'sublevelSamples', nInSub, ...
             'predictedDrop', Jhat(1) - Jhat(end));

if verbose
    fprintf('  [Phi_predictor] Upsilon=%.4g, rho=%.4f, L_J=%.4g (from %d in-sublevel samples), J_LQR(K_d)=%.4g\n', ...
        Upsilon, rho, L_J, nInSub, J_LQR_Kd);
    fprintf('  [Phi_predictor] h* = %.3f;  h_t in [%d, %d] (max graph distance %d)\n', ...
        h_star, min(h_traj), max(h_traj), maxD);
    fprintf('  [Phi_predictor] certified generations (h_{t-1} > h*): %d / %d\n', ...
        out.nCertifiedGens, maxGen);
    if out.nCertifiedGens == 0
        fprintf('  [Phi_predictor] *** Theorem 1''s hypothesis never holds: Theorem 1 certifies ZERO decrease.\n');
        fprintf('  [Phi_predictor]     L_J would have to be < %.4g (it is %.4g) for h* < max h_t.\n', ...
            L_J_crit, L_J);
    end
    fprintf('  [Phi_predictor] Phi(h_0)=%.4g vs w_c=%.4g, max P_imp=%.4g, predicted total drop=%.4g\n', ...
        Phi_traj(1), w_c, max(P_imp), out.predictedDrop);
end

end

%% ======================== Helpers ========================
function E = accumarray_blocks(K_s, colNode, numNodes)
% Collapse the Nu x Nx support pattern into an Nu x numNodes block-presence
% matrix: entry (u, j) counts nonzeros of row u in the columns owned by node j.
E = zeros(size(K_s, 1), numNodes);
for j = 1:numNodes
    cols = (colNode == j);
    if any(cols)
        E(:, j) = sum(K_s(:, cols), 2);
    end
end
end

function val = get_opt(s, field, default)
if isfield(s, field)
    val = s.(field);
else
    val = default;
end
end
