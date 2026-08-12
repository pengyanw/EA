function val = get_lqr_cost(A, B, Q, R, K, Sigma)
%GET_LQR_COST  Infinite-horizon LQR cost of the closed loop A+BK.
%
%   J(K) = E_{x_0 ~ N(0,Sigma)} [ sum_k x_k'Q x_k + u_k'R u_k ]
%        = trace(P*Sigma),
%   where P solves the discrete Lyapunov equation
%        (A+BK)' P (A+BK) - P + Q + K'RK = 0.
%   Returns Inf when A+BK is not Schur stable.
%
%   Sigma defaults to the identity, i.e. J(K) = trace(P).
%
%   NOTE: this previously returned norm(P,'fro')^2 = sum(lambda_i(P).^2)
%   rather than trace(P) = sum(lambda_i(P)). Both are increasing in P, so the
%   EA's ranking was broadly similar, but they are different objectives: the
%   Frobenius version weights the largest closed-loop mode quadratically and so
%   grows far faster near the stability boundary. That inflated the Lipschitz
%   constant L_J used in Section IV by roughly one order of magnitude
%   (L_J/J_LQR(K_d): 11.3 -> 7.5 on the 5x5 grid, 8.2 -> 5.0 on the 7x7).

    if nargin < 6, Sigma = eye(size(A,1)); end

    % Reject marginally stable closed loops. A bare rho < 1 test is not enough:
    % on plants with rho(A) = 1 exactly (e.g. the IEEE 13-bus) the EA reaches
    % gains with 1 - rho ~ 1e-15, which pass the test but make the Lyapunov
    % solve numerically meaningless -- dlyap then returns a non-PSD P with
    % trace(P) ~ -4e15. Under the old ||P||_F^2 objective that surfaced as a huge
    % POSITIVE cost and was discarded; under trace(P) it is a huge NEGATIVE cost,
    % so elitism locks onto it and the whole run is destroyed. The true cost
    % diverges as rho -> 1 anyway, so Inf is the right answer here.
    rhoTol = 1e-8;
    if any(abs(eig(A + B*K)) >= 1 - rhoTol)
        val = inf;
        return;
    end

    try
        P = dlyap((A + B*K)', Q + K'*R*K);
    catch
        val = inf;
        return;
    end

    val = trace(P * Sigma);

    % Backstop: P must be PSD for Q, R >= 0, so a non-positive trace means the
    % solve failed regardless of the margin above.
    if ~isfinite(val) || val <= 0
        val = inf;
    end
end