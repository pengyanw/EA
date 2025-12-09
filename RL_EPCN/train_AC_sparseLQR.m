clear; clc;
addpath(genpath(pwd));  
close all
%% 1. System Generation (No changes here)
% =========================================================================
gridSize       = 5;
connectThresh  = 0.5;
Ts             = 0.2;
actDensity     = 1;
seed           = 5; 

numNodes    = gridSize*gridSize;
[adjMtx, nodeCoords, susceptMtx, inertiasInv, dampings] = generate_grid_topology(gridSize, connectThresh, seed);
plot_graph(adjMtx, nodeCoords, 'k');

numActs       = round(actDensity*numNodes);
actuatedNodes = randsample(numNodes, numActs);    
sys           = generate_grid_plant(actuatedNodes, adjMtx, susceptMtx, inertiasInv, dampings, Ts);

Nx = sys.Nx;
Nu = sys.Nu;
A  = sys.A;
lamda = 0.1;
%A = A + lamda*eye(size(A));
B = sys.B2;
%B= B + 1e-1*randn(size(B));
% === RL agent init (do this once near the top of script) ===
Nu = size(B, 2);
Nx = size(A, 1);
state_dim = 5;        % 与 extract_state_features 保持一致
rng(42);              % 固定随机种子（可改）
 
%% compute the graph properties
G = graph(adjMtx);
D = distances(G);
D(isinf(D)) = max(D(~isinf(D))) + 1;

n_s  = Nx / Nu;
D_for_K = kron(D, ones(1, n_s));
use_gaussian = true;
sigma  = 1.5;
beta   = 0.7;
w_min  = 0.0;

modelFile = 'cache/rf_Kimportance_simple.mat';
K_map = predict_K_importance(A, B, adjMtx, modelFile);  % RF 模型预测出的kernel
imagesc(K_map); colorbar; title('Predicted K importance');
W = K_map;  % 直接用RF输出
%% 3. Benchmark Calculation (No changes here)
% =========================================================================
Q_bm  = eye(Nx);
R_bm_ = eye(Nu);
K_bm  = -dlqr(A, B, Q_bm, R_bm_); % Benchmark dense LQR controller
costBM = get_lqr_cost(A, B, Q_bm, R_bm_, K_bm);
fprintf('Benchmark LQR cost (dense controller): %f\n', costBM);

% Benchmark: semi-truncated K
KSuppBM     = abs(K_bm) > 1e-2;
K1          = zeros(size(K_bm));
K1(KSuppBM) = K_bm(KSuppBM);
alpha = 0;
cost_bm     = cost_EA(A, B, Q_bm, R_bm_, K1, costBM, alpha)
fprintf('Benchmark cost (dense controller): %f\n', cost_bm);
K_dense = K1;
%% 

    % ==== 0. 超参数 ====
    maxIter     = 20000;
    batchSize   = 256;
    dAct        = 64;      % actor latent 维度
    stepMag     = 0.05;    % 每次最大增量
    bufSize     = 40000;
    lr_actor    = 1e-3;
    lr_critic   = 1e-3;

    % ==== 1. 系统 embedding（固定）====
    phi_sys = build_system_embedding(A, B, adjMtx);

    % ==== 2. 可编辑位置（不用全K，减少动作维度）====
    editableIdx = get_editable_indices(K_dense, adjMtx);
    disp(editableIdx)
    % ==== 3. 初始化 K ====
    K = K_dense;     % 或者用 semi-sparse K_bm
    
    % ==== 4. 初始化 Actor / Critic 参数 ====
    sDim = length(embed_state(A,B,adjMtx,K,phi_sys,K_dense));
    actor = init_actor(sDim, dAct);
    critic = init_critic(sDim, dAct);

    % Adam 优化器
    actorOpt = struct();   % 不需要 m/v!
    criticOpt = struct();

    % ==== 5. buffer ====
    buf = init_buffer(bufSize, sDim, dAct);
%% 

   for it = 1:maxIter

    
    s = embed_state(A, B, adjMtx, K, phi_sys, K_dense);

    
    z = actor_forward(actor, s);               % dAct × 1
    z = z + 0.1 * randn(size(z));              % exploration noise

    
    K_prop = apply_action_to_K(K, z, editableIdx);

   
    K_new = safe_project_to_stable(A, B, K_dense, K_prop);

    
    cost = cost_EA(A, B, Q_bm, R_bm_, K_new, costBM, alpha);
    reward = -cost;

    
    add_transition(buf, s, z, reward);

    
    if buf.count >= batchSize
        [S_batch, Z_batch, R_batch] = sample_batch(buf, batchSize);
        % S_batch: sDim × B
        % Z_batch: dAct × B
        % R_batch: 1 × B

       
        Q_pred = critic_forward(critic, S_batch, Z_batch);    % 1 × B
        target = -R_batch;                                    % Q ≈ cost
        diffQ = Q_pred - target;                              % 1 × B

        gradQ = critic_backward(critic, S_batch, Z_batch, diffQ);
        critic = adam_update(critic, gradQ, criticOpt, lr_critic);

        
        Z_actor = actor_forward(actor, S_batch);              % dAct × B
        % Critic 评价 actor 输出的动作
        Q_actor = critic_forward(critic, S_batch, Z_actor);   % 1 × B

        gradA = actor_backward(actor, critic, S_batch, Z_actor);
        actor = adam_update(actor, gradA, actorOpt, lr_actor);
    end

    
    oldCost = cost_EA(A, B, Q_bm, R_bm_, K, costBM, alpha);
    if cost < 2*oldCost
        K = K_new;
    end

    
    if mod(it, 500) == 0
        fprintf("Iter %d | cost=%.3f | sparsity=%.3f\n", ...
            it, cost, nnz(K)/numel(K));
    end

end



    

    K_new == K
    final_cost = cost_EA(A, B, Q_bm, R_bm_, K_new, costBM, alpha)
