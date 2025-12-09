function s = extract_state_features(A, B_, K, iGen, maxGen, unstable_count, popSize)
% 抽取RL状态特征：维度=5（需与 init_rl_agent 的 state_dim 一致）
% 1) ρ(A+B*K)   2) 稀疏度   3) K的能量密度   4) 稳定比例   5) 进度

    Acl = A + B_ * K;
    rho = max(abs(eig(Acl)));

    sparsity = nnz(K) / numel(K);
    energy   = norm(K, 'fro') / max(1, numel(K));
    stable_ratio = 1 - unstable_count / max(1, popSize);
    progress = iGen / max(1, maxGen);

    % 归一/裁剪，增强数值稳定
    rho_clip = min(rho, 2.0);  % 超过2就截断
    s = [rho_clip; sparsity; energy; stable_ratio; progress];
end
