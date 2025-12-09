function [mask_r, mask_c, agent] = rl_refine_cts( ...
    agent, A, B_, K_in, cts, iGen, maxGen, unstable_count, popSize)
% 用Actor-Critic产生行/列mask，并在线更新参数
% - 输入：当前K_in（Step C 后），以及统计量
% - 输出：mask_r (Nu x 1), mask_c (Nx x 1)，并更新 agent

    % === 1) 状态提取（当前） ===
    s     = extract_state_features(A, B_, K_in, iGen, maxGen, unstable_count, popSize);
    s_col = s(:);  % 列向量

    Nu = agent.Nu; Nx = agent.Nx;

    % === 2) Actor前向：连续动作 → sigmoid → mask ===
    a = agent.actorW * s_col;          % (Nu+Nx) x 1
    a = max(min(a, 10), -10);          % 防止数值爆
    a_sig = 1 ./ (1 + exp(-a));        % sigmoid

    mask_r = a_sig(1:Nu);
    mask_c = a_sig(Nu+1:end);

    % 最小强度下限，避免整行/列被硬删
    % eps_floor = agent.clamp_min;
    % mask_r = eps_floor + (1 - eps_floor) * mask_r;
    % mask_c = eps_floor + (1 - eps_floor) * mask_c;

    % === 3) 应用动作，得到新的 K ===
    K_out = K_in .* (mask_r * mask_c');

    % === 4) 计算奖励（即时）===
    % 奖励兼顾稳定性、稀疏性与控制能量（无需等外层cost）
    Acl = A + B_ * K_out;
    rho = max(abs(eig(Acl)));
    sparsity = nnz(K_out) / numel(K_out);
    energy   = norm(K_out, 'fro') / max(1, numel(K_out));

    % 你也可以把 energy 换成与基线的相对变化，或引入外层 cost_EA 的值
    r = -energy ...
        - agent.lambda_stab * max(0, rho - 1)^2 ...
        - agent.lambda_sparse * sparsity;

    % === 5) 下一个状态（用 K_out 估计）===
    s_next = extract_state_features(A, B_, K_out, iGen, maxGen, unstable_count, popSize);
    s_next_col = s_next(:);

    % === 6) Critic 更新（TD误差）===
    V_s     = agent.criticW * s_col;
    V_next  = agent.criticW * s_next_col;
    delta   = r + agent.gamma * V_next - V_s;

    agent.criticW = agent.criticW + agent.lr_critic * delta * s_col';

    % === 7) Actor 更新（用TD误差近似 advantage 的 policy gradient）===
    % 线性策略：∂a/∂W = s_col'，用 δ 当作 advantage
    agent.actorW  = agent.actorW  + agent.lr_actor  * delta * (s_col') ;

    % === 8) 输出 mask（供外层乘回K）===
    % (已在前面赋值好)
end
