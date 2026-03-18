function [result, history] = ea_lqr_codesign(A, B, ea_params, options)
% EA_LQR_CODESIGN  Evolutionary algorithm for sparse LQR co-design
%
% Syntax:
%   [result, history] = ea_lqr_codesign(A, B, ea_params)
%   [result, history] = ea_lqr_codesign(A, B, ea_params, options)
%
% Inputs:
%   A          - State matrix (Nx x Nx)
%   B          - Input matrix (Nx x Nu)
%   ea_params  - Struct with EA parameters:
%                .popSize     - Population size (default: 20)
%                .maxGen      - Maximum generations (default: 150)
%                .pMutate     - Mutation probability (default: 0.05)
%                .pCross      - Crossover probability (default: 0.8)
%                .nTop        - Number of elite individuals (default: 10)
%                .alpha       - Performance vs hardware tradeoff (default: 0)
%   options    - Struct with optional parameters:
%                .max_Q_val   - Max Q diagonal value (default: 10)
%                .max_R_val   - Max R diagonal value (default: 10)
%                .min_Q_val   - Min Q diagonal value (default: 1e-4)
%                .min_R_val   - Min R diagonal value (default: 1e-4)
%                .verbose     - Display progress (default: true)
%                .Q_bm        - Benchmark Q matrix (default: eye(Nx))
%                .R_bm        - Benchmark R matrix (default: eye(Nu))
%                .useGersRepair - Enable Gershgorin stability repair (default: true)
%
% Outputs:
%   result     - Struct containing:
%                .K_sparse    - Best sparse controller
%                .Q_best      - Best Q matrix
%                .R_best      - Best R matrix
%                .cost        - Best cost value
%                .gene        - Best gene
%                .active_actuators - Active actuator indices
%                .active_sensors   - Active sensor indices
%   history    - Struct containing:
%                .bestCost    - Best cost per generation
%                .avgCost     - Average cost per generation
%                .unstableCount - Unstable individuals per generation
%                .deltaK_norm - K perturbation norm per generation
%                .repairCount - Gershgorin repairs per generation (if enabled)
%
% Example:
%   ea_params.popSize = 20;
%   ea_params.maxGen = 100;
%   [result, history] = ea_lqr_codesign(A, B, ea_params);
%
% Author: Pengyang Wu
% Date: 2026-02-03

%% Parse inputs
[Nx, Nu] = size(B);

% Default EA parameters
if ~isfield(ea_params, 'popSize'),  ea_params.popSize = 20;    end
if ~isfield(ea_params, 'maxGen'),   ea_params.maxGen = 150;    end
if ~isfield(ea_params, 'pMutate'),  ea_params.pMutate = 0.05;  end
if ~isfield(ea_params, 'pCross'),   ea_params.pCross = 0.8;    end
if ~isfield(ea_params, 'nTop'),     ea_params.nTop = 10;       end
if ~isfield(ea_params, 'alpha'),    ea_params.alpha = 0;       end

% Default options
if nargin < 4, options = struct(); end
if ~isfield(options, 'max_Q_val'),  options.max_Q_val = 10;    end
if ~isfield(options, 'max_R_val'),  options.max_R_val = 10;    end
if ~isfield(options, 'min_Q_val'),  options.min_Q_val = 1e-4;  end
if ~isfield(options, 'min_R_val'),  options.min_R_val = 1e-4;  end
if ~isfield(options, 'verbose'),    options.verbose = true;    end
if ~isfield(options, 'Q_bm'),       options.Q_bm = eye(Nx);    end
if ~isfield(options, 'R_bm'),       options.R_bm = eye(Nu);    end
if ~isfield(options, 'useGersRepair'), options.useGersRepair = False; end

% Gene structure
len_diag_Q    = Nx;
len_diag_R    = Nu;
len_num_links = 1;
len_act_bin   = Nu;
len_sens_bin  = Nx;
geneLength    = len_diag_Q + len_diag_R + len_num_links + len_act_bin + len_sens_bin;

% Link bounds
min_links = 1;
max_links = Nu * Nx;

%% Benchmark calculation
K_bm = -dlqr(A, B, options.Q_bm, options.R_bm);
costBM = get_lqr_cost(A, B, options.Q_bm, options.R_bm, K_bm);
K_bm(abs(K_bm)<=1e-3)=0;
if options.verbose
    fprintf('========== EA-LQR Co-Design ==========\n');
    fprintf('System: Nx=%d, Nu=%d\n', Nx, Nu);
    fprintf('Benchmark LQR cost: %.6f\n', costBM);
    fprintf('EA Parameters: popSize=%d, maxGen=%d\n', ea_params.popSize, ea_params.maxGen);
    fprintf('======================================\n\n');
end

%% Initialize population
if options.verbose, fprintf('Initializing population...\n'); end

pop = cell(ea_params.popSize, 1);
for i = 1:ea_params.popSize
    diag_Q = (options.max_Q_val - options.min_Q_val) * rand(1, len_diag_Q) + options.min_Q_val;
    diag_R = (options.max_R_val - options.min_R_val) * rand(1, len_diag_R) + options.min_R_val;
    num_links = randi([min_links, max_links]);
    act_bin  = double(rand(1, len_act_bin)  > 0.01);
    sens_bin = double(rand(1, len_sens_bin) > 0.01);
    
    pop{i} = [diag_Q, diag_R, num_links, act_bin, sens_bin];
end

% History tracking
historyBestCost = zeros(ea_params.maxGen, 1);
historyAvgCost  = zeros(ea_params.maxGen, 1);
errbuffer       = zeros(ea_params.maxGen, 1);
delta_K_norm    = zeros(ea_params.maxGen, 1);
repairBuffer    = zeros(ea_params.maxGen, 1);

% Gershgorin repair options (passed to gersgorin_stabilize_K)
gers_opts.maxIter   = 80;
gers_opts.targetRho = 0.95;
gers_opts.verbose   = false;

%% Main evolutionary loop
if options.verbose, fprintf('Starting evolution...\n'); end

for iGen = 1:ea_params.maxGen
    costs = zeros(ea_params.popSize, 1);
    unstable_count = 0;
    repair_count   = 0;
    
    % --- Fitness Evaluation ---
    for i = 1:ea_params.popSize
        gene = pop{i};
        
        % Decode gene
        diag_Q    = gene(1:len_diag_Q);
        diag_R    = gene(len_diag_Q + 1 : len_diag_Q + len_diag_R);
        num_links = round(gene(len_diag_Q + len_diag_R + 1));
        act_bin   = logical(gene(len_diag_Q + len_diag_R + len_num_links + (1:len_act_bin)));
        sens_bin  = logical(gene(end - len_sens_bin + 1:end));
        
        Q = diag(diag_Q);
        R = diag(diag_R);
        
        % Compute dense controller
        try
            K_dense = -dlqr(A, B, Q, R);
        catch
            costs(i) = 1e12;
            unstable_count = unstable_count + 1;
            continue;
        end
        
        % Sparsify by strongest elements
        K_sparse = zeros(Nu, Nx);
        [~, sorted_idx] = sort(abs(K_dense(:)), 'descend');
        keep_idx = sorted_idx(1:min(num_links, end));
        K_sparse(keep_idx) = K_dense(keep_idx);
        
        % Apply masks
        K_sparse(~act_bin, :) = 0;
        K_sparse(:, ~sens_bin) = 0;
        
        % Gershgorin stability repair (if enabled)
        if options.useGersRepair
            rho_pre = max(abs(eig(A + B * K_sparse)));
            if rho_pre >= 1.0
                [K_sparse, gers_info] = gersgorin_stabilize_K( ...
                    A, B, K_sparse, K_dense, gers_opts);
                if gers_info.success
                    repair_count = repair_count + 1;
                end
            end
        end
        
        % Compute cost
        J = cost_EA(A, B, options.Q_bm, options.R_bm, K_sparse, costBM, ea_params.alpha, max_links);
        costs(i) = J;
        
        if J >= 1e9
            unstable_count = unstable_count + 1;
        end
    end
    
    errbuffer(iGen)    = unstable_count;
    repairBuffer(iGen) = repair_count;
    
    % --- Selection ---
    [sortedCosts, sortedIdx] = sort(costs);
    pop = pop(sortedIdx);
    
    bestCost = sortedCosts(1);
    historyBestCost(iGen) = bestCost;
    historyAvgCost(iGen)  = mean(costs(costs < 1e9));
    
    if options.verbose
        if options.useGersRepair
            fprintf('Gen %3d: Best=%.4f, Avg=%.4f, Unstable=%d/%d, Repaired=%d\n', ...
                iGen, bestCost, historyAvgCost(iGen), unstable_count, ea_params.popSize, repair_count);
        else
            fprintf('Gen %3d: Best=%.4f, Avg=%.4f, Unstable=%d/%d\n', ...
                iGen, bestCost, historyAvgCost(iGen), unstable_count, ea_params.popSize);
        end
    end
    
    % Compute delta K norm
    bestGene = pop{1};
    diag_Q_best = bestGene(1:len_diag_Q);
    diag_R_best = bestGene(len_diag_Q + 1 : len_diag_Q + len_diag_R);
    Q_best = diag(diag_Q_best);
    R_best = diag(diag_R_best);
    best_K_dense = -dlqr(A, B, Q_best, R_best);
    
    [~, sorted_idx] = sort(abs(best_K_dense(:)), 'descend');
    keep_idx = sorted_idx(1:min(bestGene(len_diag_Q + len_diag_R + 1), end));
    K_sparse_best = zeros(Nu, Nx);
    K_sparse_best(keep_idx) = best_K_dense(keep_idx);
    K_sparse_best(~logical(bestGene(len_diag_Q + len_diag_R + len_num_links + (1:len_act_bin))), :) = 0;
    K_sparse_best(:, ~logical(bestGene(end - len_sens_bin + 1:end))) = 0;
    
    delta_K_norm(iGen) = norm(K_sparse_best - best_K_dense, 2) / (norm(best_K_dense, 2) + eps);
    
    % --- Create next generation ---
    newPop = cell(ea_params.popSize, 1);
    newPop(1:ea_params.nTop) = pop(1:ea_params.nTop); % Elitism
    
    for i = (ea_params.nTop+1):ea_params.popSize
        % Parent selection
        parent1 = pop{randi(ea_params.popSize)};
        parent2 = pop{randi(ea_params.popSize)};
        
        % Crossover
        if rand < ea_params.pCross
            crossPoint = randi(geneLength - 1);
            child = [parent1(1:crossPoint), parent2(crossPoint+1:end)];
        else
            child = parent1;
        end
        
        % Mutation (continuous parts)
        for j = 1:(len_diag_Q + len_diag_R + len_num_links)
            if rand < ea_params.pMutate
                child(j) = child(j) + randn * 0.1 * (options.max_Q_val - options.min_Q_val);
            end
        end
        
        % Bound enforcement
        child(1:len_diag_Q) = min(max(child(1:len_diag_Q), options.min_Q_val), options.max_Q_val);
        child(len_diag_Q+1 : len_diag_Q+len_diag_R) = ...
            min(max(child(len_diag_Q+1 : len_diag_Q+len_diag_R), options.min_R_val), options.max_R_val);
        child(len_diag_Q+len_diag_R+1) = ...
            min(max(round(child(len_diag_Q+len_diag_R+1)), min_links), max_links);
        
        % Mutation (binary parts)
        for j = (len_diag_Q + len_diag_R + len_num_links + 1):geneLength
            if rand < ea_params.pMutate
                child(j) = 1 - child(j);
            end
        end
        
        newPop{i} = child;
    end
    
    pop = newPop;
end

%% Extract results
bestGene = pop{1};
diag_Q_best = bestGene(1:len_diag_Q);
diag_R_best = bestGene(len_diag_Q + 1 : len_diag_Q + len_diag_R);
num_links_best = round(bestGene(len_diag_Q + len_diag_R + 1));
act_bin_best = logical(bestGene(len_diag_Q + len_diag_R + len_num_links + (1:len_act_bin)));
sens_bin_best = logical(bestGene(end - len_sens_bin + 1:end));

Q_best = diag(diag_Q_best);
R_best = diag(diag_R_best);
K_dense_best = -dlqr(A, B, Q_best, R_best);

[~, sorted_idx] = sort(abs(K_dense_best(:)), 'descend');
keep_idx = sorted_idx(1:min(num_links_best, end));
K_sparse_best = zeros(Nu, Nx);
K_sparse_best(keep_idx) = K_dense_best(keep_idx);
K_sparse_best(~act_bin_best, :) = 0;
K_sparse_best(:, ~sens_bin_best) = 0;

% Apply Gershgorin repair to final result
if options.useGersRepair
    rho_final = max(abs(eig(A + B * K_sparse_best)));
    if rho_final >= 1.0
        [K_sparse_best, ~] = gersgorin_stabilize_K( ...
            A, B, K_sparse_best, K_dense_best, gers_opts);
    end
end

% Pack results
result.K_sparse = K_sparse_best;
result.Q_best = Q_best;
result.R_best = R_best;
result.cost = historyBestCost(end);
result.gene = bestGene;
result.active_actuators = find(act_bin_best);
result.active_sensors = find(sens_bin_best);
result.num_links = nnz(K_sparse_best);
cost_BM = get_lqr_cost(A, B, options.Q_bm, options.R_bm, K_bm);
result.costBM = cost_EA(A, B, options.Q_bm, options.R_bm, K_bm, cost_BM, ea_params.alpha, max_links);
result.improvement = (result.costBM - result.cost) / result.costBM;

history.bestCost = historyBestCost;
history.avgCost = historyAvgCost;
history.unstableCount = errbuffer;
history.deltaK_norm = delta_K_norm;
history.repairCount = repairBuffer;
history.generations = 1:ea_params.maxGen;

%% Summary
if options.verbose
    fprintf('\n========== EA Summary ==========\n');
    fprintf('Best cost: %.6f\n', result.cost);
    fprintf('Benchmark cost: %.6f\n', result.costBM);
    fprintf('Improvement: %.2f%%\n', result.improvement * 100);
    fprintf('Active actuators: %d/%d\n', length(result.active_actuators), Nu);
    fprintf('Active sensors: %d/%d\n', length(result.active_sensors), Nx);
    fprintf('Controller links: %d/%d\n', result.num_links, Nu*Nx);
    fprintf('Max delta K/K: %.4f\n', max(delta_K_norm));
    fprintf('================================\n');
end

end
