function agent = init_rl_agent(Nu, Nx, state_dim)
% 轻量 Actor-Critic（线性）初始化
% action 维度 = Nu(行mask) + Nx(列mask)

    act_dim = Nu + Nx;

    agent.Nu = Nu;
    agent.Nx = Nx;
    agent.state_dim = state_dim;
    agent.act_dim   = act_dim;

    % 线性 actor / critic 权重
    agent.actorW  = 0.01 * randn(act_dim, state_dim);  % a = W*s
    agent.criticW = 0.01 * randn(1,       state_dim);  % V = w*s

    % 超参数（可调）
    agent.gamma         = 0.95;   % 折扣
    agent.lr_actor      = 0.01;   % 学习率（actor）
    agent.lr_critic     = 0.01;   % 学习率（critic）
    agent.clamp_min     = 0.20;   % mask 不低于20%强度，防止整行/整列骤删
    agent.lambda_stab   = 50;     % 稳定性罚系数
    agent.lambda_sparse = 0.10;   % 稀疏性罚系数
end
