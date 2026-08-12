function [gapTraj, out] = gap_predictor(A, B, K_dense, geneHist, w, opts)
%GAP_PREDICTOR  Certified optimality gap of Theorem C / Corollary C.1 (proof_IV.md).
%
%   At each generation the elite gene theta*_t = [ell, a, s] determines a retained
%   actuator set S* (the support of a). Holding ell and s fixed, write
%
%       m_u  = links actuator u contributes  = nnz of row u of Pi_ell(K_d)|_s
%       M(T) = w_a|T| + w_c * sum_{u in T} m_u              (exactly modular)
%
%   Theorem C: if the LQR value f(S) = -J_LQR(K(S))/J_LQR(K_d) is monotone on the
%   feasible family and has submodularity ratio gamma, then for any competitor S,
%
%       U(S) - U(S*) <= M(S* \ S) + (1/gamma - 1) * M(S \ S*).
%
%   S_opt is unknown, but S*\S_opt is a subset of S* and S_opt\S* of Omega\S*, so
%
%       min_S J_EA(S) >= J_EA(S*) - [ M(S*) + (1/gamma - 1) * M(Omega\S*) ]
%
%   and the bracket -- returned here as gapTraj -- is computable from the run.
%
%   IMPORTANT: m_u must be evaluated at the CURRENT generation's ell and s. Reusing
%   N_c from generation 1 (a different sensor mask) inflates M(Omega\S*) by ~20x.
%
%   Scope (see proof_IV.md section 10): the bound is stated for the ACTUATOR
%   coordinate at fixed ell and s. The sensor value function is -inf on most
%   subsets, so no analogous bound is available there; Theorem A (a.s. convergence
%   to a 1-flip local optimum) does cover both coordinates.
%
%   SLACK TERM. Theorem C's proof uses the ADDITION half of 1-flip local
%   optimality, f(S*+u) - f(S*) <= w_a + w_c*m_u for all u not in S*. The elite
%   does not satisfy this at every generation -- it can prune past a local optimum
%   mid-run, and on the 7x7 grid it has not reached one by generation 150. Rather
%   than gate the bound on that check, we carry the violation explicitly:
%
%       eps_u := max(0, f(S*+u) - f(S*) - w_a - w_c*m_u)     (profit of adding u)
%
%   Repeating the proof with f(S*+u) - f(S*) <= w_a + w_c*m_u + eps_u gives
%
%       min_a J_EA >= J_EA(S*) - [ M(S*) + (1/g - 1) M(Omega\S*) + (1/g) sum_u eps_u ]
%
%   which holds at EVERY generation and collapses to the original expression
%   exactly when the elite is 1-flip locally optimal (all eps_u = 0). The eps_u
%   come from the same Lyapunov solves the local-optimality check needs, so there
%   is no extra cost, and sum(eps_u) measures how far the run is from converging.
%
% Syntax:
%   [gapTraj, out] = gap_predictor(A, B, K_dense, geneHist, w)
%   [gapTraj, out] = gap_predictor(..., opts)
%
% Inputs:
%   A, B      - plant matrices
%   K_dense   - K_d, already thresholded (the matrix the EA sorts and truncates)
%   geneHist  - maxGen x (1+Nu+Nx) elite genes, from history.bestGene
%   w         - struct with fields wc, wa (ws is not used: sensors are held fixed)
%   opts      - optional:
%               .gamma     submodularity ratio; if absent it is estimated here
%               .nSample   DR triples for the gamma estimate (default 250)
%               .pct       percentile of the DR ratio taken as gamma (default 5)
%               .Q, .R     LQR weights (default identity)
%               .verbose   print the constants (default true)
%
% Outputs:
%   gapTraj - maxGen x 1 certified gap in raw J_EA units
%   out     - struct: gamma, drViol, M_S, M_C, Na, Nc, gammaSamples,
%             valid (maxGen x 1 logical), nViol (maxGen x 1 count of profitable
%             single-actuator additions), validFrac, validAtEnd
%
% Author: Pengyang Wu

if nargin < 6, opts = struct(); end
Q       = getf(opts, 'Q',       eye(size(A,1)));
R       = getf(opts, 'R',       eye(size(B,2)));
nSample = getf(opts, 'nSample', 250);
pct     = getf(opts, 'pct',     5);
verbose = getf(opts, 'verbose', true);

[Nu, Nx] = size(K_dense);
maxGen   = size(geneHist, 1);
idx_a    = 1 + (1:Nu);
idx_s    = 1 + Nu + (1:Nx);

[~, sortIdx] = sort(abs(K_dense(:)), 'descend');
nnzKd = nnz(K_dense);

%% ---------- submodularity ratio of f ----------
if isfield(opts, 'gamma')
    gamma  = opts.gamma;
    drViol = NaN;
    ratios = [];
else
    % Sample gamma where Theorem C USES it. The theorem is invoked at S = S*,
    % the returned architecture, whose size is 8-30% of N_u; sampling the
    % 45-80% band (the previous default) measures a region neither search ever
    % occupies, and measures it badly: there the aggregate gain f(S u U) - f(S)
    % is near zero, so the quotient is dominated by its denominator. Measured
    % per band on the 5x5 and 7x7 grids, min over samples:
    %     |S|/N_u    0.05-0.15   0.15-0.30   0.30-0.45   0.45-0.80
    %     gamma_min      >0.976      -0.119      -2.012       0.110
    % One gamma is used for the whole curve, so the early-generation bound is
    % approximate; it is dominated by M(S*) there in any case.
    % Take the band from the TAIL of the trajectory, not all of it. The search
    % starts at the all-ones mask, so spanning the whole run would put the dense
    % end -- where gamma is worst -- back into the estimate, which is exactly the
    % region the restriction is meant to exclude. Theorem C is invoked at
    % S = S*, so the last quarter of the run is the relevant range.
    tail = max(1, floor(0.75*maxGen)) : maxGen;
    kTail = zeros(numel(tail),1);
    for tt = 1:numel(tail)
        kTail(tt) = nnz(logical(geneHist(tail(tt), idx_a)));
    end
    band = [max(2, floor(0.75*min(kTail))), min(Nu-2, ceil(1.25*max(kTail)))];
    if band(2) <= band(1)
        band = [max(2, min(kTail)), min(Nu-2, max(kTail)+2)];
    end
    if band(2) <= band(1), band = [max(2,Nu-3), Nu-2]; end   % degenerate fallback
    [gamma, drViol, ratios] = estimate_gamma(A, B, Q, R, K_dense, sortIdx, ...
                                             nnzKd, Nu, Nx, nSample, pct, band);
end
if ~(gamma > 0)
    error('gap_predictor:gamma', ...
        'Estimated gamma = %.4g is not positive; Theorem C does not apply.', gamma);
end

%% ---------- per-generation gap and validity ----------
checkValid = getf(opts, 'checkValid', true);
gapTraj = zeros(maxGen, 1);
M_S = zeros(maxGen, 1);  M_C = zeros(maxGen, 1);
Na  = zeros(maxGen, 1);  Nc  = zeros(maxGen, 1);
valid = true(maxGen, 1);  nViol = zeros(maxGen, 1);  epsSum = zeros(maxGen, 1);
J_Kd = get_lqr_cost(A, B, Q, R, K_dense);

for t = 1:maxGen
    ell = max(1, min(round(geneHist(t, 1)), nnzKd));
    av  = logical(geneHist(t, idx_a));
    sv  = logical(geneHist(t, idx_s));

    % Pi_ell(K_d) masked to the CURRENT sensor set, all actuators available
    Kpool = zeros(Nu, Nx);
    Kpool(sortIdx(1:ell)) = K_dense(sortIdx(1:ell));
    Kpool(:, ~sv) = 0;
    m = sum(Kpool ~= 0, 2);              % m_u for every u in Omega

    % M must be charged on the DECODED controller, not on the mask. An actuator
    % whose row is emptied by the ell-truncation costs nothing in J_EA -- one does
    % not install an actuator with no controller entries -- so only rows with
    % m_u > 0 contribute the w_a term. The two agree whenever ell is large enough
    % that every masked row survives, which is why the EA (ell ~ 90-99% of max)
    % never exposed the difference and the greedy baseline (ell ~ 5-12%) does: at
    % one greedy point 16 actuators were masked on but only 4 had any entries,
    % inflating M(S*) from 2.0 to 6.8 and the gap from 2.06 to 6.82.
    %
    % The per-element weight w_a*1{m_u>0} + w_c*m_u is still fixed given (ell, s),
    % so M stays exactly modular and Lemma 1 / Theorem C' are unaffected.
    liveS  = av  & (m' > 0);
    liveC  = ~av & (m' > 0);
    M_S(t) = w.wa * nnz(liveS)  + w.wc * sum(m(av));
    M_C(t) = w.wa * nnz(liveC)  + w.wc * sum(m(~av));
    Na(t)  = nnz(liveS);
    Nc(t)  = sum(m(av));

    if checkValid
        % Addition half of 1-flip local optimality, carried as a slack term rather
        % than used as a gate. An addition that destabilises satisfies the
        % inequality trivially (f(S*+u) = -inf), so it contributes eps_u = 0.
        Kc = Kpool;  Kc(~av, :) = 0;
        Jc = get_lqr_cost(A, B, Q, R, Kc);
        if ~isfinite(Jc)
            valid(t) = false;  nViol(t) = NaN;  epsSum(t) = Inf;
            gapTraj(t) = Inf;                   % elite unstable: nothing certified
            continue;
        end
        fS = -Jc / J_Kd;  v = 0;  es = 0;
        for u = find(~av)
            K2 = Kc;  K2(u, :) = Kpool(u, :);
            J2 = get_lqr_cost(A, B, Q, R, K2);
            if ~isfinite(J2), continue; end
            e = (-J2/J_Kd) - fS - w.wa - w.wc * m(u);
            if e > 1e-9, v = v + 1;  es = es + e; end
        end
        nViol(t) = v;  valid(t) = (v == 0);  epsSum(t) = es;
    end

    gapTraj(t) = M_S(t) + (1/gamma - 1) * M_C(t) + (1/gamma) * epsSum(t);
end

out = struct('gamma', gamma, 'drViolations', drViol, 'M_S', M_S, 'M_C', M_C, ...
             'Na', Na, 'Nc', Nc, 'gammaSamples', ratios, ...
             'valid', valid, 'nViol', nViol, 'epsSum', epsSum, ...
             'validFrac', mean(valid), 'validAtEnd', valid(end));

if verbose
    fprintf('  [gap_predictor] gamma = %.4f (DR violations %.1f%%), 1/gamma-1 = %.4f\n', ...
        gamma, 100*drViol, 1/gamma - 1);
    fprintf('  [gap_predictor] terminal: N_a=%d, N_c=%d | M(S*)=%.3f + (1/g-1)M(Om\\S*)=%.3f + (1/g)sum(eps)=%.3f = gap %.3f\n', ...
        Na(end), Nc(end), M_S(end), (1/gamma-1)*M_C(end), epsSum(end)/gamma, gapTraj(end));
    if checkValid
        fprintf('  [gap_predictor] elite is 1-flip locally optimal in %d/%d generations; at the end: %s (%d profitable additions, slack %.3f)\n', ...
            nnz(valid), maxGen, string(valid(end)), nViol(end), epsSum(end));
    end
end
end

%% ======================== helpers ========================
function [gamma, viol, ratios] = estimate_gamma(A, B, Q, R, Kd, sortIdx, nz, Nu, Nx, nSample, pct, band)
% Das-Kempe submodularity ratio of f(S) = -J_LQR(K(S))/J_LQR(K_d):
%
%     gamma = min over (S, U)  sum_{u in U} [f(S+u) - f(S)]  /  [f(S∪U) - f(S)]
%
% NOT the pairwise diminishing-returns ratio [f(S+v)-f(S)]/[f(T+v)-f(T)]. That
% proxy has a SINGLE marginal in the denominator, so it explodes (and changes sign)
% whenever that marginal is near zero: it returned 0.986 on one seed and -17 on
% another for the same plant family. The Das-Kempe denominator is an aggregate gain
% and is far better conditioned -- across five seeds on each grid, the per-sample
% ratio has median exactly 1.0000 and per-seed minimum 0.88-0.99.
%
% Sampled on the band the search actually occupies (passed in), with |U| >= 2.
% Returns a PERCENTILE (default 5th) rather than the sample minimum; on the
% operating band the two nearly coincide (min >= 0.976 across the tested
% instances), unlike on the 45-80% band where the minimum falls to 0.110.
J_Kd = get_lqr_cost(A, B, Q, R, Kd);
Kfull = zeros(Nu, Nx);
Kfull(sortIdx(1:nz)) = Kd(sortIdx(1:nz));
f = @(S) fval(A, B, Q, R, Kfull, J_Kd, S, Nu, Nx);

rng(7);                                  % fixed: the estimate must be reproducible
ratios = []; valid = 0; viol = 0;
for k = 1:nSample
    p  = randperm(Nu);
    kS = randi(band);                                  % the band the search occupies
    % |U| >= 2. With |U| = 1 the numerator has a single term equal to the
    % denominator, so the quotient is identically 1 and carries no information
    % about submodularity; including such draws inflates the estimate.
    kU = randi([2, max(3, round(0.15*Nu))]);
    S    = p(1:kS);
    rest = setdiff(1:Nu, S);
    if numel(rest) < kU, continue; end
    U = rest(randperm(numel(rest), kU));

    fS = f(S);  fSU = f(union(S, U));
    if ~isfinite(fS) || ~isfinite(fSU), continue; end
    den = fSU - fS;
    if den <= 1e-8, continue; end                     % need a genuine aggregate gain

    num = 0; ok = true;
    for u = U
        fu = f(union(S, u));
        if ~isfinite(fu), ok = false; break; end
        num = num + (fu - fS);
    end
    if ~ok, continue; end

    valid = valid + 1;
    ratios(end+1) = num / den;                                     %#ok<AGROW>
    if num < den - 1e-9, viol = viol + 1; end
end
if valid < 20
    error('gap_predictor:gammaSamples', ...
        'Only %d valid Das-Kempe samples; gamma cannot be estimated here.', valid);
end
gamma = min(1, prctile(ratios, pct));
viol  = viol / valid;
end

function y = fval(A, B, Q, R, Kfull, J_Kd, S, Nu, Nx)
K = zeros(Nu, Nx);
K(S, :) = Kfull(S, :);
Jl = get_lqr_cost(A, B, Q, R, K);
if ~isfinite(Jl), y = -inf; else, y = -Jl / J_Kd; end
end

function v = getf(s, f, d)
if isfield(s, f), v = s.(f); else, v = d; end
end
