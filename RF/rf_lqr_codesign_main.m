
    clear; clc; close all;
    addpath(genpath(pwd));

    %% 1) System generation (same)
    gridSize       = 3;
    connectThresh  = 0.5;
    Ts             = 0.2;
    actDensity     = 1;
    seed           = 8;

    numNodes    = gridSize*gridSize;
    [adjMtx, nodeCoords, susceptMtx, inertiasInv, dampings] = ...
        generate_grid_topology(gridSize, connectThresh, seed);
    figure; plot_graph(adjMtx, nodeCoords, 'k'); title('Grid topology');

    numActs       = round(actDensity*numNodes);
    actuatedNodes = randsample(numNodes, numActs);
    sys           = generate_grid_plant(actuatedNodes, adjMtx, susceptMtx, inertiasInv, dampings, Ts);

    Nx = sys.Nx;   Nu = sys.Nu;
    A  = sys.A;    B_ = sys.B2;

    %% 2) Base LQR and cost
    Q_bm  = eye(Nx);
    R_bm_ = eye(Nu);
    K_bm  = -dlqr(A, B_, Q_bm, R_bm_);
    costBM = get_lqr_cost(A, B_, Q_bm, R_bm_, K_bm);
    max_links = Nu * Nx;

    
    alpha_cost = 0.5;            
    stab_thresh = 1.0;            % stability mark

    fprintf('Benchmark LQR cost: %.6f\n', costBM);
    cost_bm = alpha_cost*1 + (1-alpha_cost)*nnz(K_bm)/max_links;
    fprintf('Benchmark blended cost: %.6f\n', cost_bm);

    %% 3) Gene structure
    len_diag_Q    = Nx;
    len_diag_R    = Nu;
    len_num_links = 1;
    geneLength    = len_diag_Q + len_diag_R + len_num_links;

    
    max_Q_val = 10; min_Q_val = 1e-4;
    max_R_val = 10; min_R_val = 1e-4;
    min_links = 1;  % at least 1
    % max_links 

    %% 4) 构建训练集（从随机基因 + LQR 稀疏化得到 K）
    Ntrain = 200;
    gene_mat = zeros(Ntrain, geneLength);
    K_cell   = cell(Ntrain,1);

    fprintf('Building training set with %d samples...\n', Ntrain);
    keep = false(Ntrain,1);
    s = 0;
    while s < Ntrain
        % 随机基因
        diag_Q    = (max_Q_val - min_Q_val)*rand(1, len_diag_Q) + min_Q_val;
        diag_R    = (max_R_val - min_R_val)*rand(1, len_diag_R) + min_R_val;
        num_links = randi([min_links, max_links]);
        gene = [diag_Q, diag_R, num_links];

        % 解码并求稠密 LQR
        Q = diag(diag_Q); R = diag(diag_R);
        try
            K_dense = -dlqr(A, B_, Q, R);
        catch
            continue; % R 非正定等失败，跳过
        end

        % 稀疏化
        K_sparse = sparsify_by_links(K_dense, num_links);

        % 稳定性检查与损失（需要稳定才纳入训练）
        [J, isStable, ~] = eval_K_cost(K_sparse, A, B_, Q_bm, R_bm_, costBM, alpha_cost, max_links);
        if ~isStable || ~isfinite(J), continue; end

        % 记录
        s = s + 1;
        gene_mat(s,:) = gene;
        K_cell{s}     = K_sparse;
        keep(s)       = true;
    end
    gene_mat = gene_mat(keep,:);
    K_cell   = K_cell(keep);
    Ntrain   = size(gene_mat,1);
    fprintf('Training set size: %d\n', Ntrain);

    %% 5) 训练 RF 模型（每个 K(i,j) 一棵回归森林）
    fprintf('Training Random Forest models...\n');
    models = rf_train_K_models(gene_mat, K_cell, Nu, Nx);

    %% 6) 主循环：RF 生成候选、评估、主动学习式迭代
    rounds     = 50;         % iterations
    batchSize  = 20;         % every round, # of trials
    historyBestCost = nan(rounds,1);
    historyAvgCost  = nan(rounds,1);
    historyUnstable = nan(rounds,1);

    bestOverallCost = inf; bestGene = []; bestK = [];

    for r = 1:rounds
        J_batch = nan(batchSize,1);
        stableMask = false(batchSize,1);

        bestJ_round = inf; bestK_round = []; bestGene_round = [];

        for b = 1:batchSize
            % 采样一个新基因（可替换为基因邻域扰动/贝叶斯优化等）
            diag_Q    = (max_Q_val - min_Q_val)*rand(1, len_diag_Q) + min_Q_val;
            diag_R    = (max_R_val - min_R_val)*rand(1, len_diag_R) + min_R_val;
            num_links = randi([min_links, max_links]);
            gene_new  = [diag_Q, diag_R, num_links];

            % RF 预测 K 并评估
            [K_new, J_new, isStable, ~] = rf_generate_once( ...
                models, gene_new, A, B_, Q_bm, R_bm_, costBM, alpha_cost, max_links);

            J_batch(b)    = J_new;
            stableMask(b) = isStable;

            % 记录轮次最优
            if isStable && J_new < bestJ_round
                bestJ_round   = J_new;
                bestK_round   = K_new;
                bestGene_round= gene_new;
            end
        end

        % 统计
        if any(stableMask)
            historyBestCost(r) = min(J_batch(stableMask));
            historyAvgCost(r)  = mean(J_batch(stableMask));
        else
            historyBestCost(r) = inf;
            historyAvgCost(r)  = inf;
        end
        historyUnstable(r) = sum(~stableMask);

        fprintf('Round %d: Best=%.4f, Avg=%.4f, Unstable=%d/%d\n', ...
            r, historyBestCost(r), historyAvgCost(r), historyUnstable(r), batchSize);

        % 主动学习：把轮次最优加入训练集并立刻重训
        if isfinite(bestJ_round)
            gene_mat = [gene_mat; bestGene_round];
            K_cell   = [K_cell;   {bestK_round}];
            models   = rf_train_K_models(gene_mat, K_cell, Nu, Nx);
        end

        % 全局最优
        if isfinite(bestJ_round) && bestJ_round < bestOverallCost
            bestOverallCost = bestJ_round;
            bestGene = bestGene_round;
            bestK    = bestK_round;
        end
    end

    %% 7) 结果与绘图（风格与现有一致）
    fprintf('\n===== RF Code-sign Summary =====\n');
    fprintf('Best blended cost: %.6f\n', bestOverallCost);
    fprintf('Advantage over benchmark: %.4f\n', (cost_bm - bestOverallCost)/cost_bm);

    if ~exist('figures','dir'), mkdir figures; end

    figure; 
    plot(1:rounds, historyBestCost, 'b-', 'LineWidth', 2); hold on;
    plot(1:rounds, historyAvgCost,  'g--', 'LineWidth', 1.5);
    grid on; xlabel('Round'); ylabel('Blended Cost (normalized)');
    title(sprintf('RF LQR Co-Design (grid=%d)', gridSize));
    legend('Best per round','Average (stable only)','Location','best');
    saveas(gcf, sprintf('figures/rf_evo_cost_grid%d.png', gridSize));

    figure;
    plot(1:rounds, historyUnstable, 'r-o', 'LineWidth', 1.5);
    grid on; xlabel('Round'); ylabel('# Unstable per batch');
    title('Unstable count per round');
    saveas(gcf, sprintf('figures/rf_unstable_grid%d.png', gridSize));

    % 可选保存最优 K
    save(sprintf('figures/rf_bestK_grid%d.mat', gridSize), 'bestK', 'bestGene', ...
         'historyBestCost','historyAvgCost','historyUnstable');

    fprintf('Saved plots and best K.\n');
