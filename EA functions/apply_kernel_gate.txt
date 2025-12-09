function [K_final, isStable] = apply_kernel_gate(K_dense, num_links, W, A, B_)
%--------------------------------------------------------------------------
%  apply_kernel_gate  --  Apply distance-based gating (Gaussian / exponential)
%  combined with top-N sparsification.
%
%  Inputs:
%     K_dense   : Dense LQR gain matrix (Nu x Nx)
%     num_links : Number of nonzero entries to keep (scalar)
%     W         : Distance-based kernel weight matrix (Nu x Nx)
%     A, B_     : System matrices
%
%  Outputs:
%     K_final   : Sparsified and kernel-gated gain
%     isStable  : Boolean flag, true if A + B_*K_final is stable
%
%--------------------------------------------------------------------------

    % 1️⃣ Apply kernel weights (Gaussian or exponential)
    K_weighted = K_dense .* W;

    % 2️⃣ Flatten and sort by absolute weighted magnitude
    [~, idx] = sort(abs(K_weighted(:)), 'descend');

    % 3️⃣ Keep top num_links entries (clipped to available elements)
    num_links = min(num_links, numel(K_dense));
    keep_idx = idx(1:num_links);

    % 4️⃣ Construct sparse K with kernel shaping
    K_final = zeros(size(K_dense));
    K_final(keep_idx) = K_weighted(keep_idx);

    % 5️⃣ Optional normalization (preserve overall scale)
    norm_dense = norm(K_dense, 'fro') + 1e-12;
    norm_final = norm(K_final, 'fro') + 1e-12;
    K_final = K_final * (norm_dense / norm_final);

    % 6️⃣ Stability check with fallback line search
    Acl = A + B_ * K_final;
    if max(abs(eig(Acl))) < 1
        isStable = true;
        return;
    end

    % --- fallback (progressively mix toward K_dense) ---
    eta = 1.0;
    while max(abs(eig(A + B_ * ((1-eta)*K_dense + eta*K_final)))) >= 1 && eta > 1e-3
        eta = eta / 2;
    end
    K_final = (1-eta)*K_dense + eta*K_final;
    isStable = max(abs(eig(A + B_*K_final))) < 1;
end
