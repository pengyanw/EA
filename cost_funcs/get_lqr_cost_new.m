function val = get_lqr_cost_new(A, B, Q, R, K)
%GET_LQR_COST  Fast surrogate of LQR Lyapunov cost (no dlyap).
%   Original code:
%       P = dlyap((A+BK)', Q+K'RK); val = ||P||_F^2
%   New surrogate:
%       estimate tr(P) via Hutchinson + finite-horizon rollout,
%       and return (tr(P))^2 as a cheap proxy for ||P||_F^2.
%
%   Complexity: O(S * T * nnz(A+BK)) instead of Lyapunov solve.
%
%   Tunable hyperparameters (set here; you can also pass via globals).

    % ----------------- hyperparameters -----------------
    T  = 80;     % rollout horizon (20~80 typical)
    S  = 10;      % Hutchinson samples (3~10 typical)
    growth_tol = 1e6;  % early divergence threshold (stability screening)
    % ---------------------------------------------------

    Acl = A + B*K;                % closed-loop
    Qbar = Q + K'*R*K;            % stage matrix in Lyapunov

    % Quick (cheap) stability screening: if rollout blows up, treat as unstable.
    % This avoids eig() and avoids dlyap() failure modes.
    tr_est = 0;
    for s = 1:S
        % Rademacher vector z in {+1,-1}^n, E[zz']=I
        z = 2*(rand(size(A,1),1) > 0.5) - 1;

        x = z;
        accum = 0;

        for t = 1:T
            % accumulate x' Qbar x  (this equals z' P_T z)
            accum = accum + (x' * Qbar * x);

            % propagate
            x = Acl * x;

            % early exit if unstable-like growth
            if ~isfinite(norm(x)) || norm(x) > growth_tol
                val = inf;
                return;
            end
        end

        tr_est = tr_est + accum;
    end

    tr_est = tr_est / S;

    % Proxy for ||P||_F^2: cheap and monotone in practice for EA ranking
    val = tr_est^2;

    % Optional: if numerical weirdness
    if ~isfinite(val)
        val = inf;
    end
end
