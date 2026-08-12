function [best, hist] = greedy_codesign(A, B, ea_params, options)
%GREEDY_CODESIGN  Round-robin greedy local search over the co-design gene.
%
%   The baseline for "why an evolutionary algorithm rather than something
%   simpler?". It searches the SAME gene theta = [ell, a, s] as
%   ea_lqr_codesign_gershgorin.m, scores it through the SAME cost path, and
%   terminates at a 1-flip local optimum of the mask coordinates -- which is
%   exactly the class of point Theorem A says the EA converges to. That makes the
%   two directly comparable rather than merely both "reasonable".
%
%   Starting from the dense architecture theta_0 = [nnz(K_d), 1, 1] (the EA's own
%   initialisation), cycle three phases, each run to its own fixed point:
%
%     1. link      : try ell +- 1 and ell +- 5, take the best improvement, repeat.
%                    The +-5 step is deliberate -- it is precisely the EA's link
%                    neighbourhood (delta ~ Unif{-d..d}, d = 5), so greedy is not
%                    handed a global scan the EA structurally cannot perform.
%     2. actuator  : flip each a_u in BOTH directions, take the best single move,
%                    repeat to a fixed point.
%     3. sensor    : same for each s_j.
%
%   Stop when a full cycle yields no improvement.
%
%   Fairness controls, all deliberate:
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
%   [best, hist] = greedy_codesign(A, B, ea_params, options)
%
% Inputs:
%   A, B       - plant matrices
%   ea_params  - struct; only .popSize (for the budget axis) and .alpha are read
%   options    - .useGersRepair (default false), .denseTol (default 1e-3),
%                .Q, .R (default identity), .maxCycles (default 20),
%                .linkSteps (default [-5 -1 1 5]), .verbose (default false)
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
linkSteps = getf(options, 'linkSteps',     [-5 -1 1 5]);
verbose   = getf(options, 'verbose',       false);
alpha     = getf(ea_params, 'alpha',       0);
Np        = getf(ea_params, 'popSize',     20);
maxGen    = getf(ea_params, 'maxGen',      150);

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

% --- initialisation: the dense architecture, as in the EA ---
ell = max_links;  a = true(1,Nu);  s = true(1,Nx);
J = oracle(ell, a, s);
record();

cycles = 0;
for c = 1:maxCycles
    cycles = c;
    J0 = J;

    % ---- phase 1: link count, within the EA's own +-d neighbourhood ----
    improved = true;
    while improved
        improved = false;
        for d = linkSteps
            e2 = ell + d;
            if e2 < 1 || e2 > max_links, continue; end
            J2 = oracle(e2, a, s);
            if J2 < J - 1e-12, J = J2; ell = e2; improved = true; end
            record();
        end
    end

    % ---- phase 2: actuators, single flips in both directions ----
    improved = true;
    while improved
        improved = false;  bestJ = J;  bu = 0;
        for u = 1:Nu
            a2 = a;  a2(u) = ~a2(u);
            J2 = oracle(ell, a2, s);
            if J2 < bestJ - 1e-12, bestJ = J2; bu = u; end
            record();
        end
        if bu > 0, a(bu) = ~a(bu); J = bestJ; improved = true; end
    end

    % ---- phase 3: sensors, single flips in both directions ----
    improved = true;
    while improved
        improved = false;  bestJ = J;  bj = 0;
        for j = 1:Nx
            s2 = s;  s2(j) = ~s2(j);
            J2 = oracle(ell, a, s2);
            if J2 < bestJ - 1e-12, bestJ = J2; bj = j; end
            record();
        end
        if bj > 0, s(bj) = ~s(bj); J = bestJ; improved = true; end
    end

    if J >= J0 - 1e-12, break; end
end

K = decode(ell, a, s);
best = struct('K_sparse', K, 'gene', [ell, double(a), double(s)], 'J', J, ...
              'Na', nnz(any(K~=0,2)), 'Ns', nnz(any(K~=0,1)), 'Nc', nnz(K), ...
              'ell', ell, 'evals', evals, 'nInf', nInf, 'cycles', cycles);

% Incumbent sampled on the EA's budget axis: generation g of the EA has consumed
% Np*g evaluations, so greedy's incumbent after Np*g evaluations is the
% like-for-like comparison. Greedy converges early, so hold its final value.
JperGen = zeros(maxGen, 1);
for g = 1:maxGen
    k = find(evalIdx <= Np*g, 1, 'last');
    if isempty(k), JperGen(g) = incumbent(1); else, JperGen(g) = incumbent(k); end
end
hist = struct('evalIdx', evalIdx, 'incumbent', incumbent, 'JperGen', JperGen);

if verbose
    fprintf('  [greedy] J=%.4f  ell=%d  N_a=%d N_s=%d N_c=%d  evals=%d (%.2f x %d-gen EA)  inf=%d  cycles=%d\n', ...
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
