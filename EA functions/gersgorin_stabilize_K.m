function [K_out, info] = gersgorin_stabilize_K(A, B, K_sparse, K_ref, opts)
% GERSGORIN_STABILIZE_K  Stabilize a sparse controller using Gershgorin disk theory.
%
%   Given a sparse K that may yield an unstable closed-loop A+B*K, this
%   function adjusts the nonzero entries to achieve Schur stability (|lambda|<1).
%
%   Three-phase approach:
%     Phase 1: Subgradient descent on max Gershgorin row-sum (preserves sparsity)
%     Phase 2: Spectral radius check; line search on sparsity support toward K_ref
%     Phase 3: Line search toward full K_ref + re-sparsification (fallback)
%
% Syntax:
%   [K_out, info] = gersgorin_stabilize_K(A, B, K_sparse, K_ref)
%   [K_out, info] = gersgorin_stabilize_K(A, B, K_sparse, K_ref, opts)
%
% Inputs:
%   A          - State matrix (Nx x Nx)
%   B          - Input matrix (Nx x Nu)
%   K_sparse   - Sparse controller gain (Nu x Nx), possibly unstable
%   K_ref      - Reference stabilizing controller (Nu x Nx), e.g., from dlqr
%   opts       - (optional) struct:
%                .maxIter    - Max subgradient iterations (default: 80)
%                .targetRho  - Target Gershgorin row-sum (default: 0.95)
%                .verbose    - Print progress (default: false)
%
% Outputs:
%   K_out      - Stabilized controller
%   info       - Diagnostics struct:
%                .phase      - Which phase succeeded (0/1/2/3, -1 if failed)
%                .success    - true if K_out is Schur stable
%                .iters      - Subgradient iterations used
%                .rho_in     - Input spectral radius
%                .rho_out    - Output spectral radius
%                .gers_in    - Input max Gershgorin row-sum
%                .gers_out   - Output max Gershgorin row-sum
%
% Theory:
%   The Gershgorin row-sum R_i = sum_j |Acl(i,j)| satisfies:
%     max_i R_i < 1  ==>  rho(Acl) < 1   (sufficient for Schur stability)
%   R_i(K) is convex in K, so the feasible set {K : max_i R_i < 1} is convex.
%   The subgradient of R_{i*} w.r.t. K(u,j) is  sign(Acl(i*,j)) * B(i*,u),
%   where i* = argmax_i R_i.  Polyak step size is used for convergence.
%
% Author: Pengyang Wu
% Date: 2026-02-03

%% Parse options
if nargin < 5, opts = struct(); end
maxIter   = get_opt(opts, 'maxIter',   80);
targetRho = get_opt(opts, 'targetRho', 0.95);
verbose   = get_opt(opts, 'verbose',   false);
% Step-size cap on the Polyak step. Inf reproduces Eq. (37) of the paper
% exactly, which is what Proposition 5 is proved for. Set to 0.5 to restore the
% previous clamped step (which is the damped Polyak step with
% lambda_t = min(1, 0.5*||g||_F^2/(Rbar - gamma)) in (0,1]).
etaCap    = get_opt(opts, 'etaCap',    Inf);

[Nu, Nx] = size(K_sparse);
mask = (K_sparse ~= 0);
nLinks = nnz(K_sparse);

%% Initial diagnostics
Acl0    = A + B * K_sparse;
rho_in  = max(abs(eig(Acl0)));
gers_in = max(sum(abs(Acl0), 2));

info = struct('phase', -1, 'success', false, 'iters', 0, ...
              'rho_in', rho_in, 'rho_out', rho_in, ...
              'gers_in', gers_in, 'gers_out', gers_in);

%% Phase 0: Already stable
if rho_in < 1.0
    K_out = K_sparse;
    info.phase   = 0;
    info.success = true;
    if verbose, fprintf('[GersStab] Already stable (rho=%.4f)\n', rho_in); end
    return;
end

%% ========== Phase 1: Subgradient descent on max Gershgorin row-sum ==========
K = K_sparse;
best_K   = K;
best_rho = rho_in;

for iter = 1:maxIter
    Acl = A + B * K;
    R   = sum(abs(Acl), 2);         % Gershgorin row-sums (Nx x 1)
    [rmax, i_star] = max(R);

    % Converged within Gershgorin bound
    if rmax < targetRho
        best_K = K;
        best_rho = max(abs(eig(A + B * K)));
        break;
    end

    % Subgradient: dR_{i*}/dK(u,j) = sign(Acl(i*,j)) * B(i*,u)

    % At Acl(i_star,j) == 0 the row-sum is non-differentiable; sign(0) = 0 is a
    % valid subgradient element since 0 in [-|B(i_star,u)|, |B(i_star,u)|].
    s_row = sign(Acl(i_star, :));    % 1 x Nx
    b_row = B(i_star, :);            % 1 x Nu
    grad  = b_row' * s_row;          % Nu x Nx (outer product)
    grad  = grad .* mask;            % preserve sparsity

    gn = norm(grad, 'fro');
    if gn < 1e-14, break; end

    % Polyak step size, Eq. (37): eta_t = (Rbar(K) - gamma) / ||g~||_F^2
    eta = (rmax - targetRho) / (gn^2 + 1e-12);
    eta = min(eta, etaCap);         % etaCap = Inf => pure Polyak step (37)

    K = K - eta * grad;

    % Proposition 5 assumes G intersect K_S is non-empty. When it is not, the
    % pure Polyak step has no descent guarantee and the iterate can diverge, so
    % stop before Inf/NaN reaches the eigenvalue solve below. Only best_K is used
    % after this loop, and it only ever holds finite iterates, so breaking is
    % enough -- no need to write K.
    if ~all(isfinite(K(:))), break; end

    % Track best by actual spectral radius (every 10 iters, cheap enough)
    if mod(iter, 10) == 0
        rho_curr = max(abs(eig(A + B * K)));
        if rho_curr < best_rho
            best_rho = rho_curr;
            best_K   = K;
        end
    end
end

info.iters = min(iter, maxIter);

% Final eigenvalue check for Phase 1 result
rho_p1 = max(abs(eig(A + B * best_K)));
if rho_p1 < best_rho
    best_rho = rho_p1;
end

if best_rho < 1.0
    K_out = best_K;
    info.phase   = 1;
    info.success = true;
    info.rho_out  = best_rho;
    info.gers_out = max(sum(abs(A + B * best_K), 2));
    if verbose
        fprintf('[GersStab] Phase 1 success: rho %.4f -> %.4f (%d iters)\n', ...
            rho_in, best_rho, info.iters);
    end
    return;
end

%% ========== Phase 2: Line search on sparsity support toward K_ref ==========
K_ref_proj = K_ref .* mask;         % project K_ref onto sparsity pattern

% Binary search: K_try = (1-eta)*best_K + eta*K_ref_proj
eta_lo = 0;
eta_hi = 1;
K_out  = K_sparse;                  % fallback
found  = false;

for bs = 1:25
    eta_mid = (eta_lo + eta_hi) / 2;
    K_try   = (1 - eta_mid) * best_K + eta_mid * K_ref_proj;
    rho_try = max(abs(eig(A + B * K_try)));

    if rho_try < 0.999
        found  = true;
        K_out  = K_try;
        eta_hi = eta_mid;           % try less K_ref, preserve more structure
    else
        eta_lo = eta_mid;
    end
end

if found
    info.phase   = 2;
    info.success = true;
    info.rho_out  = max(abs(eig(A + B * K_out)));
    info.gers_out = max(sum(abs(A + B * K_out), 2));
    if verbose
        fprintf('[GersStab] Phase 2 success: rho %.4f -> %.4f (eta=%.3f)\n', ...
            rho_in, info.rho_out, eta_hi);
    end
    return;
end

%% ========== Phase 3: Line search toward full K_ref + re-sparsification ==========
eta_lo = 0;
eta_hi = 1;
found  = false;

for bs = 1:25
    eta_mid = (eta_lo + eta_hi) / 2;
    K_full  = (1 - eta_mid) * best_K + eta_mid * K_ref;

    % Re-sparsify to original link count
    K_sp = zeros(Nu, Nx);
    [~, idx] = sort(abs(K_full(:)), 'descend');
    keep = idx(1:min(nLinks, end));
    K_sp(keep) = K_full(keep);

    rho_try = max(abs(eig(A + B * K_sp)));

    if rho_try < 0.999
        found  = true;
        K_out  = K_sp;
        eta_hi = eta_mid;
    else
        eta_lo = eta_mid;
    end
end

if found
    info.phase   = 3;
    info.success = true;
    info.rho_out  = max(abs(eig(A + B * K_out)));
    info.gers_out = max(sum(abs(A + B * K_out), 2));
    if verbose
        fprintf('[GersStab] Phase 3 success: rho %.4f -> %.4f (eta=%.3f)\n', ...
            rho_in, info.rho_out, eta_hi);
    end
    return;
end

%% All phases failed
K_out = K_sparse;
info.phase   = -1;
info.success = false;
info.rho_out  = rho_in;
info.gers_out = gers_in;
if verbose
    fprintf('[GersStab] All phases failed. rho=%.4f unchanged.\n', rho_in);
end

end

%% ======================== Helper ========================
function val = get_opt(s, field, default)
    if isfield(s, field)
        val = s.(field);
    else
        val = default;
    end
end
