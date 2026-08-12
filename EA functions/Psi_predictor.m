function [Jhat, out] = Psi_predictor(A, B, K_dense, geneHist, J0, w, ea_params, opts)
%PSI_PREDICTOR  Per-generation cost prediction from the actuator/sensor removal
%   certificate, the mask-coordinate counterpart of Phi_predictor.
%
%   At gene theta = [ell, a, s] with K = K_s(theta), removing actuator row u is
%   certified to decrease J_EA when
%
%       Delta_a(u) := w_a + w_c*|supp(K(u,:))|                (exact saving)
%             >
%       Psi_a(u)   := (L_J/J_LQR(K_d))*( sigma_crit*||K(u,:)||_F
%                                        + 0.5*||K(u,:)||_F^2 )
%
%   which follows from the descent lemma with ||grad J(K)||_F <= L_J*||K-K_d||_F
%   <= L_J*sigma_crit (Theorem 2's stability condition) and grad J(K_d) = 0.
%   Sensors are identical with w_s and column norms.
%
%   With a dedicated elite-child slot that mutates exactly one uniformly chosen
%   mask bit, the probability of landing on a certified, currently-active bit is
%   c_t/(N_u+N_x), so
%
%       E[ J_EA(t-1) - J_EA(t) ] >= (c_t/(N_u+N_x)) * min_{v in cert}(Delta-Psi),
%
%   and the recursion below accumulates that bound along the MEASURED gene
%   trajectory (the bound is a per-generation drift given the current state; it
%   does not predict the trajectory itself).
%
% Inputs:
%   A, B      - plant matrices
%   K_dense   - K_d (already thresholded; same matrix the EA truncates)
%   geneHist  - maxGen x (1+Nu+Nx) elite genes, from history.bestGene
%   J0        - measured J_EA of the elite at generation 1 (raw units)
%   w         - struct with fields wc, wa, ws
%   ea_params - struct with popSize, nTop
%   opts      - optional: .Q .R .verbose .sublevel .statesPerNode
%
% Outputs:
%   Jhat - predicted cost trajectory (upper envelope of the measured curve)
%   out  - struct with sigma_crit, L_J, Upsilon, rho, cTraj, gainTraj, driftTraj
%
% Author: Pengyang Wu

if nargin < 8, opts = struct(); end
Q       = getf(opts, 'Q', eye(size(A,1)));
R       = getf(opts, 'R', eye(size(B,2)));
verbose = getf(opts, 'verbose', true);

[Nu, Nx] = size(K_dense);
n        = Nu + Nx;
maxGen   = size(geneHist, 1);
idx_nL   = 1;  idx_a = 1+(1:Nu);  idx_s = 1+Nu+(1:Nx);

%% ---- constants ----
J_LQR_Kd = get_lqr_cost(A, B, Q, R, K_dense);
[~, sortIdx] = sort(abs(K_dense(:)), 'descend');
nnzKd = nnz(K_dense);

% Upsilon, rho from Lemma 1 (upper envelope of |K_d,ij| vs graph distance)
[Ups, rho] = decayFit(K_dense, opts, Nu, Nx);

% sigma_crit, Eq. (23). L from Assumption 1(a).
L = max([norm(A), norm(B), norm(Q), norm(R), 1 + 1e-9]);
sigma_crit = (1 - rho^2) / (4 * Ups^2 * L * (L + 2*Ups*rho));

% L_J on the pruning family, restricted to a sublevel set (Lemma 2)
subMult = getf(opts, 'sublevel', 5);
L_J = 0;
for ell = unique(round(linspace(1, nnzKd, 120)))
    Ks = zeros(Nu, Nx);  Ks(sortIdx(1:ell)) = K_dense(sortIdx(1:ell));
    d2 = norm(Ks - K_dense, 'fro')^2;
    if d2 < eps, continue; end
    Jl = get_lqr_cost(A, B, Q, R, Ks);
    if ~isfinite(Jl) || Jl > subMult * J_LQR_Kd, continue; end
    L_J = max(L_J, 2 * (Jl - J_LQR_Kd) / d2);
end
kappa = L_J / J_LQR_Kd;

%% ---- per-generation certificate ----
cTraj     = zeros(maxGen, 1);
gainTraj  = zeros(maxGen, 1);
cAddTraj  = zeros(maxGen, 1);
gSumTraj  = zeros(maxGen, 1);
for t = 1:maxGen
    ell = max(1, min(round(geneHist(t, idx_nL)), nnzKd));
    K = zeros(Nu, Nx);  K(sortIdx(1:ell)) = K_dense(sortIdx(1:ell));
    K(~logical(geneHist(t, idx_a)), :) = 0;
    K(:, ~logical(geneHist(t, idx_s))) = 0;

    rN = vecnorm(K, 2, 2);   rM = sum(K ~= 0, 2);
    cN = vecnorm(K, 2, 1)';  cM = sum(K ~= 0, 1)';

    Da = w.wa + w.wc * rM;   Pa = kappa * (sigma_crit * rN + 0.5 * rN.^2);
    Ds = w.ws + w.wc * cM;   Ps = kappa * (sigma_crit * cN + 0.5 * cN.^2);

    okA = (rM > 0) & (Da > Pa);      % certified AND currently active
    okS = (cM > 0) & (Ds > Ps);

    cTraj(t) = nnz(okA) + nnz(okS);
    g = [Da(okA) - Pa(okA); Ds(okS) - Ps(okS)];
    if isempty(g), gainTraj(t) = 0; else, gainTraj(t) = min(g); end

    % --- additive gains, for the multi-flip expectation bound ---
    % Rows have pairwise disjoint supports, and so do columns; only a row and a
    % column can share an entry. Claiming the link saving for rows only makes
    % the per-bit gains exactly additive over any subset:
    %     g_add(u) = w_a + w_c|row u| - Psi_a(u)   (actuators)
    %     g_add(j) = w_s             - Psi_s(j)    (sensors)
    % Psi is sub-additive since ||Delta K_S||_F^2 <= sum_v ||v||^2 and
    % sqrt(sum x^2) <= sum x, so Psi(S) <= sum_v Psi(v).
    gaA = Da - Pa;                      % rows already exclude double counting
    gaS = (w.ws - Ps);                  % drop the link term for columns
    aA  = (rM > 0) & (gaA > 0);
    aS  = (cM > 0) & (gaS > 0);
    cAddTraj(t)  = nnz(aA) + nnz(aS);
    gSumTraj(t)  = sum(gaA(aA)) + sum(gaS(aS));
end

%% ---- recursion ----
% Bound A (single uniformly chosen bit; needs a dedicated single-bit elite slot)
driftA = (cTraj / n) .* gainTraj;

% Bound B (the paper's per-bit Bernoulli operator, many bits may flip at once).
% Condition on every bit outside the strictly certified set staying put:
%     E[dJ] >= (1-p_m)^(n-c') * p_m * sum_{v in C'} g_add(v)
pm     = getf(ea_params, 'pMutate', 0.05);
driftB = ((1-pm).^(n - cAddTraj)) .* pm .* gSumTraj;

% Both are valid lower bounds, so take the better one each generation.
driftTraj = max(driftA, driftB);
Jhat = zeros(maxGen, 1);  Jhat(1) = J0;
for t = 2:maxGen
    Jhat(t) = max(Jhat(t-1) - driftTraj(t-1), 0);
end

out = struct('sigma_crit', sigma_crit, 'L_J', L_J, 'kappa', kappa, ...
             'Upsilon', Ups, 'rho', rho, 'J_LQR_Kd', J_LQR_Kd, ...
             'cTraj', cTraj, 'gainTraj', gainTraj, ...
             'cAddTraj', cAddTraj, 'gSumTraj', gSumTraj, ...
             'driftA', driftA, 'driftB', driftB, 'driftTraj', driftTraj, ...
             'predictedDrop', Jhat(1) - Jhat(end), ...
             'dropA', sum(driftA(1:end-1)), 'dropB', sum(driftB(1:end-1)), ...
             'nCertifiedGens', nnz(cTraj > 0), 'winB', nnz(driftB > driftA));

if verbose
    fprintf('  [Psi_predictor] Upsilon=%.4g rho=%.4f L_J=%.4g L_J/J=%.4g sigma_crit=%.4g\n', ...
        Ups, rho, L_J, kappa, sigma_crit);
    fprintf('  [Psi_predictor] certified bits c_t: %d -> %d (of %d);  generations with c_t>0: %d/%d\n', ...
        cTraj(1), cTraj(end), n, out.nCertifiedGens, maxGen);
    fprintf('  [Psi_predictor] single-bit bound  A: total %.4f\n', out.dropA);
    fprintf('  [Psi_predictor] multi-flip bound  B: total %.4f   (B wins in %d/%d gens)\n', ...
        out.dropB, out.winB, maxGen);
    fprintf('  [Psi_predictor] max(A,B) per-gen drift: %.4f -> %.4f;  total %.4f\n', ...
        driftTraj(1), driftTraj(end), out.predictedDrop);
end
end

%% ======================== helpers ========================
function [Ups, rho] = decayFit(Kd, opts, Nu, Nx)
    spn = getf(opts, 'statesPerNode', 2);
    adj = getf(opts, 'adjMtx', []);
    if isempty(adj)
        error('Psi_predictor:adj', 'opts.adjMtx is required to fit Lemma 1.');
    end
    D = distances(graph(adj ~= 0));
    % State-to-node map. The grid plants interleave (node ii owns [2ii-1, 2ii]);
    % the IEEE 13-bus blocks (node ii owns [ii, N+ii]). Pass opts.colNode for
    % anything that is not interleaved.
    colNode = getf(opts, 'colNode', ceil((1:Nx)/spn));
    if numel(colNode) ~= Nx
        error('Psi_predictor:colNode', ...
            'opts.colNode must have %d entries, got %d.', Nx, numel(colNode));
    end
    Dk = D((1:Nu)', colNode);
    hV = unique(Dk(isfinite(Dk)))';  mM = nan(size(hV));
    for k = 1:numel(hV)
        v = abs(Kd(Dk == hV(k)));  v = v(v > 0);
        if ~isempty(v), mM(k) = max(v); end
    end
    ok = isfinite(mM) & mM > 0;
    cf = polyfit(hV(ok), log(mM(ok)), 1);
    rho = min(max(exp(cf(1)), 1e-6), 1 - 1e-6);
    Ups = max(1, exp(cf(2)));
    Ups = Ups * max(1, max(abs(Kd(:)) ./ (Ups * rho.^Dk(:))));
end

function v = getf(s, f, d)
    if isfield(s, f), v = s.(f); else, v = d; end
end
