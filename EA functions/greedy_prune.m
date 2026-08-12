function [best, hist] = greedy_prune(A, B, ea_params, options)
%GREEDY_PRUNE  Monotone (pruning-only) greedy over the co-design gene.
%
%   The B2 baseline: "why an evolutionary algorithm rather than something
%   simpler?". Same gene theta = [ell, a, s] as ea_lqr_codesign_gershgorin.m,
%   same decode, same cost oracle -- but the move set is restricted to
%   REMOVALS only. Nothing that has been pruned is ever restored:
%
%     1. link      : ell <- ell - delta, delta in {1, 5} only. Never increases.
%                    The -5 step is the EA's own link neighbourhood
%                    (delta ~ Unif{-d..d}, d = 5), so greedy is not handed a
%                    global scan the EA structurally cannot perform.
%     2. actuator  : a_u : 1 -> 0 only, best single removal, to a fixed point.
%     3. sensor    : s_j : 1 -> 0 only, best single removal, to a fixed point.
%
%   Cycle the three phases; stop when a full cycle yields no improvement.
%   Starting from the dense architecture theta_0 = [nnz(K_d), 1, 1] (the EA's
%   own initialisation), the retained support is therefore nested decreasing
%   over the whole run -- this is classical backward elimination, and it is the
%   procedure the "greedy pruning" discussion in Section IV actually refers to.
%
%   Difference from GREEDY_CODESIGN (kept for reference): that variant flips
%   masks in BOTH directions and allows ell to increase, so it can backtrack and
%   terminates at a 1-flip local optimum. GREEDY_PRUNE cannot backtrack, so its
%   fixed point is only "no single removal improves" -- a strictly weaker
%   stopping condition. Consequently the slack Xi in Theorem 4 need not vanish
%   at a greedy_prune output, whereas it can at a greedy_codesign output.
%
%   Fairness controls, all deliberate and identical to GREEDY_CODESIGN:
%     - one cost oracle, transcribed from ea_lqr_codesign_gershgorin.m:206-243,
%       including the Gershgorin repair branch when options.useGersRepair is set;
%     - NO memoisation. The EA re-evaluates its elites every generation, so
%       caching greedy's repeated probes would be an asymmetric discount;
%     - destabilising moves are SCORED (as 1e9), not skipped, so the budget
%       reflects what the search actually costs;
%     - every oracle call is counted, and the incumbent is recorded against that
%       count, so greedy and the EA can be compared on a common budget axis
%       (the EA spends exactly N_p evaluations per generation).
%
% Syntax:
%   [best, hist] = greedy_prune(A, B, ea_params, options)
%
% Inputs:
%   A, B       - plant matrices
%   ea_params  - struct; only .popSize, .maxGen (budget axis) and .alpha are read
%   options    - .useGersRepair (default false), .denseTol (default 1e-3),
%                .Q, .R (default identity), .maxCycles (default 20),
%                .linkSteps (default [5 1], magnitudes; all applied as
%                decrements), .verbose (default false)
%
% Outputs:
%   best - struct: K_sparse, gene [ell a s], J, Na, Ns, Nc, evals, nInf, cycles
%   hist - struct: evalIdx (n x 1), incumbent (n x 1), and JperGen (maxGen x 1),
%          the incumbent sampled at e = popSize * g so it can be plotted directly
%          against the EA's generation axis
%
% Author: Pengyang Wu

if nargin < 4, options = struct(); end
Nx = size(A,1);  Nu = size(B,2);
Q         = getf(options, 'Q',             eye(Nx));
R         = getf(options, 'R',             eye(Nu));
denseTol  = getf(options, 'denseTol',      1e-3);
useRepair = getf(options, 'useGersRepair', false);
maxCycles = getf(options, 'maxCycles',     20);
verbose   = getf(options, 'verbose',       false);
alpha     = getf(ea_params, 'alpha',       0);
Np        = getf(ea_params, 'popSize',     20);
maxGen    = getf(ea_params, 'maxGen',      150);

% The link step is the EA's own mutation range d, read from the SAME field the
% EA reads (ea_lqr_codesign_gershgorin.m: ea_params.mutRange, default 5). Do not
% hard-code it here: Algorithm 3's input line claims "the same d as Algorithm 1"
% and that claim has to survive someone changing d.
d         = getf(ea_params, 'mutRange',    5);
linkSteps = getf(options, 'linkSteps',     [d 1]);
skipNoops = getf(options, 'skipNoops',     true);

% Pruning-only: link moves are decrements whatever sign the caller passed.
linkSteps = -abs(linkSteps(:)');

w_c = 0.05*(1-alpha);  w_a = 0.40*(1-alpha);  w_s = 0.20*(1-alpha);

% --- EA conventions, copied so the two searches decode identically ---
K_dense = -dlqr(A, B, Q, R);
K_dense(abs(K_dense) <= denseTol) = 0;
costBM  = get_lqr_cost(A, B, Q, R, K_dense);
[~, K_sorted_idx] = sort(abs(K_dense(:)), 'descend');
max_links = nnz(K_dense);
gers_opts = struct('etaCap', Inf);

evals = 0;  nInf = 0;
evalIdx = zeros(0,1);  incumbent = zeros(0,1);

% Committed moves and the probes spent to commit them, per phase [link a s].
% The three phases differ sharply on both counts, so keeping them apart is the
% only way to read the flat steps in the Fig. 1 greedy curve.
nMove = zeros(1,3);  nProbe = zeros(1,3);

% Probes spent on a mask bit that is on but whose row/column is already empty in
% the decoded controller. Clearing such a bit cannot change K, so the probe is a
% provable no-op -- it is charged a full Lyapunov solve and can never be
% selected. Diagnostic only; the search is not altered.
nNoop = zeros(1,3);

% --- initialisation: the dense architecture, as in the EA ---
ell = max_links;  a = true(1,Nu);  s = true(1,Nx);
J = oracle(ell, a, s);
record();

cycles = 0;
for c = 1:maxCycles
    cycles = c;
    J0 = J;

    % ---- phase 1: link count, decrements only ----
    improved = true;
    while improved
        improved = false;  bestJ = J;  bell = ell;
        for d = linkSteps
            e2 = ell + d;
            if e2 < 1, continue; end
            J2 = oracle(e2, a, s);
            nProbe(1) = nProbe(1) + 1;
            if J2 < bestJ - 1e-12, bestJ = J2; bell = e2; end
            record();
        end
        if bell < ell, ell = bell; J = bestJ; improved = true; nMove(1) = nMove(1)+1; end
    end

    % ---- phase 2: actuators, 1 -> 0 only ----
    improved = true;
    while improved
        improved = false;  bestJ = J;  bu = 0;
        Kcur = decode(ell, a, s);
        for u = find(a)
            % Clearing a mask bit whose row is already empty cannot change K, so
            % the probe is a provable no-op: it can never be selected, and
            % detecting it needs no Lyapunov solve. Skipping leaves the sequence
            % of committed moves identical and only lowers the evaluation count.
            % Symmetric to the EA reusing its elites' costs; charging one and not
            % the other would bias the shared budget axis.
            if skipNoops && ~any(Kcur(u,:)), nNoop(2) = nNoop(2) + 1; continue; end
            a2 = a;  a2(u) = false;
            J2 = oracle(ell, a2, s);
            nProbe(2) = nProbe(2) + 1;
            if ~any(Kcur(u,:)), nNoop(2) = nNoop(2) + 1; end
            if J2 < bestJ - 1e-12, bestJ = J2; bu = u; end
            record();
        end
        if bu > 0, a(bu) = false; J = bestJ; improved = true; nMove(2) = nMove(2)+1; end
    end

    % ---- phase 3: sensors, 1 -> 0 only ----
    improved = true;
    while improved
        improved = false;  bestJ = J;  bj = 0;
        Kcur = decode(ell, a, s);
        for j = find(s)
            if skipNoops && ~any(Kcur(:,j)), nNoop(3) = nNoop(3) + 1; continue; end
            s2 = s;  s2(j) = false;
            J2 = oracle(ell, a, s2);
            nProbe(3) = nProbe(3) + 1;
            if ~any(Kcur(:,j)), nNoop(3) = nNoop(3) + 1; end
            if J2 < bestJ - 1e-12, bestJ = J2; bj = j; end
            record();
        end
        if bj > 0, s(bj) = false; J = bestJ; improved = true; nMove(3) = nMove(3)+1; end
    end

    if J >= J0 - 1e-12, break; end
end

K = decode(ell, a, s);
best = struct('K_sparse', K, 'gene', [ell, double(a), double(s)], 'J', J, ...
              'Na', nnz(any(K~=0,2)), 'Ns', nnz(any(K~=0,1)), 'Nc', nnz(K), ...
              'ell', ell, 'evals', evals, 'nInf', nInf, 'cycles', cycles);

% Incumbent sampled on the EA's budget axis. The caller supplies the EA's
% cumulative evaluation count per generation; without it we assume the EA pays
% N_p per generation. When the EA reuses its elites' costs it pays N_p in the
% first generation and only N_p - n_e afterwards, which roughly halves the axis
% and must be reflected here or the two curves are not on a common budget.
eaEvals = getf(options, 'eaEvalSchedule', Np * (1:maxGen)');
JperGen = zeros(maxGen, 1);
for g = 1:maxGen
    k = find(evalIdx <= eaEvals(g), 1, 'last');
    if isempty(k), JperGen(g) = incumbent(1); else, JperGen(g) = incumbent(k); end
end
hist = struct('evalIdx', evalIdx, 'incumbent', incumbent, 'JperGen', JperGen, ...
              'nMove', nMove, 'nProbe', nProbe, 'nNoop', nNoop);

if verbose
    fprintf('  [prune] J=%.4f  ell=%d  N_a=%d N_s=%d N_c=%d  evals=%d (%.2f x %d-gen EA)  inf=%d  cycles=%d\n', ...
        best.J, best.ell, best.Na, best.Ns, best.Nc, evals, evals/(Np*maxGen), maxGen, nInf, cycles);
end

%% ===================== nested helpers =====================
    function K = decode(e, av, sv)
        e = max(1, min(round(e), max_links));
        K = zeros(Nu, Nx);
        keep = K_sorted_idx(1:e);
        K(keep) = K_dense(keep);
        K(~av, :) = 0;
        K(:, ~sv) = 0;
    end

    function Jv = oracle(e, av, sv)
        evals = evals + 1;
        K = decode(e, av, sv);
        if useRepair && max(abs(eig(A + B*K))) >= 1.0
            [K, ~] = gersgorin_stabilize_K(A, B, K, K_dense, gers_opts);
        end
        Jl = get_lqr_cost(A, B, Q, R, K);
        if ~isfinite(Jl)
            Jv = 1e9;  nInf = nInf + 1;  return;
        end
        Jv = Jl/costBM + w_a*nnz(any(K~=0,2)) + w_s*nnz(any(K~=0,1)) + w_c*nnz(K);
    end

    function record()
        evalIdx(end+1,1)   = evals;   %#ok<AGROW>
        incumbent(end+1,1) = J;       %#ok<AGROW>
    end
end

function v = getf(s, f, d)
if isstruct(s) && isfield(s, f), v = s.(f); else, v = d; end
end
