function [result, history] = ea_lqr_codesign_gershgorin(A, B, ea_params, options)
% EA_LQR_CODESIGN_NEW  Evolutionary sparse LQR co-design (improved)
%
% Key differences from ea_lqr_codesign:
%   1) Q, R fixed to benchmark values (dlqr computed once, not per individual)
%   2) Gene = [n_links, b_act, b_sens] -- no Q/R in chromosome
%   3) Softmax parent selection: P(i) = exp(-J_i / tau) / Z
%   4) Segment-aware crossover:
%        Integer  (n_links) : fitness-weighted interpolation
%        Binary   (b_act, b_sens) : fitness-weighted per-bit uniform
%   5) 3-phase Gershgorin hard repair during evolution
%      (gersgorin_stabilize_K applied to every unstable individual)
%
% Syntax:
%   [result, history] = ea_lqr_codesign_new(A, B, ea_params)
%   [result, history] = ea_lqr_codesign_new(A, B, ea_params, options)
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
%                .verbose       - Display progress (default: true)
%                .Q_bm          - Fixed Q for LQR (default: eye(Nx))
%                .R_bm          - Fixed R for LQR (default: eye(Nu))
%                .gers_weight   - Gershgorin soft penalty weight in cost_EA
%                                 (default: 0). Set >0 to add soft penalty.
%                .useGersRepair - 3-phase hard repair on every unstable
%                                 individual during evolution (default: true)
%
% Outputs:
%   result     - Struct containing:
%                .K_sparse         - Best sparse controller
%                .Q_best           - Q matrix used (= Q_bm)
%                .R_best           - R matrix used (= R_bm)
%                .cost             - Best cost value
%                .gene             - Best gene [n_links, b_act, b_sens]
%                .active_actuators - Active actuator indices
%                .active_sensors   - Active sensor indices
%                .num_links        - Number of non-zero elements in K
%                .costBM           - Benchmark cost
%                .improvement      - Relative improvement over benchmark
%   history    - Struct containing:
%                .bestCost      - Best cost per generation
%                .avgCost       - Average cost per generation
%                .unstableCount - Unstable individuals per generation
%                .deltaK_norm   - K perturbation norm per generation
%                .generations   - Generation indices
%
% Example:
%   ea_params.popSize = 20;
%   ea_params.maxGen  = 100;
%   [result, history] = ea_lqr_codesign_new(A, B, ea_params);
%
% Author: Pengyang Wu
% Date: 2026-02-10

%% ==================== Parse inputs ====================
[Nx, Nu] = size(B);

% EA parameters
if ~isfield(ea_params, 'popSize'),  ea_params.popSize = 20;    end
if ~isfield(ea_params, 'maxGen'),   ea_params.maxGen  = 150;   end
if ~isfield(ea_params, 'pMutate'),  ea_params.pMutate = 0.05;  end
if ~isfield(ea_params, 'pCross'),   ea_params.pCross  = 0.8;   end
if ~isfield(ea_params, 'nTop'),     ea_params.nTop    = 10;    end
if ~isfield(ea_params, 'alpha'),    ea_params.alpha   = 0;     end

% Options
if nargin < 4, options = struct(); end
if ~isfield(options, 'verbose'),       options.verbose       = true;  end
if ~isfield(options, 'Q_bm'),          options.Q_bm          = eye(Nx); end
if ~isfield(options, 'R_bm'),          options.R_bm          = eye(Nu); end
if ~isfield(options, 'gers_weight'),   options.gers_weight   = 0;     end
if ~isfield(options, 'useGersRepair'), options.useGersRepair = true;   end

% Gene structure: [n_links(1), b_act(Nu), b_sens(Nx)]
idx_nL   = 1;
idx_act  = 1 + (1:Nu);
idx_sens = 1 + Nu + (1:Nx);
geneLength = 1 + Nu + Nx;

% Link bounds
min_links = 1;
max_links = Nu * Nx;

%% ==================== Benchmark (computed once) ====================
Q_fixed = options.Q_bm;
R_fixed = options.R_bm;
K_dense = -dlqr(A, B, Q_fixed, R_fixed);


K_bm = K_dense;
K_bm(abs(K_bm) <= 1e-3) = 0;
costBM  = get_lqr_cost(A, B, Q_fixed, R_fixed, K_bm);
% Pre-sort K_dense magnitude (shared across all individuals)
[~, K_sorted_idx] = sort(abs(K_dense(:)), 'descend');

if options.verbose
    fprintf('========== EA-LQR Co-Design (New) ==========\n');
    fprintf('System: Nx=%d, Nu=%d\n', Nx, Nu);
    fprintf('Q, R fixed to benchmark (dlqr computed once)\n');
    fprintf('Benchmark LQR cost: %.6f\n', costBM);
    fprintf('Gene: [n_links, b_act(%d), b_sens(%d)], length=%d\n', Nu, Nx, geneLength);
    fprintf('EA: popSize=%d, maxGen=%d, pMut=%.3f, pCross=%.2f\n', ...
        ea_params.popSize, ea_params.maxGen, ea_params.pMutate, ea_params.pCross);
    fprintf('Gershgorin hard repair: %s (soft penalty weight: %.1f)\n', ...
        mat2str(options.useGersRepair), options.gers_weight);
    fprintf('Crossover: softmax selection + segment-aware\n');
    fprintf('============================================\n\n');
end

%% ==================== Initialize population ====================
if options.verbose, fprintf('Initializing population...\n'); end

pop = cell(ea_params.popSize, 1);
for i = 1:ea_params.popSize
    num_links = randi([min_links, max_links]);
    act_bin   = double(rand(1, Nu) > 0.01);
    sens_bin  = double(rand(1, Nx) > 0.01);
    pop{i} = [num_links, act_bin, sens_bin];
end

% History tracking
historyBestCost = zeros(ea_params.maxGen, 1);
historyAvgCost  = zeros(ea_params.maxGen, 1);
errbuffer       = zeros(ea_params.maxGen, 1);
repairBuffer    = zeros(ea_params.maxGen, 1);
delta_K_norm    = zeros(ea_params.maxGen, 1);

% Gershgorin hard repair options
gers_opts.maxIter   = 80;
gers_opts.targetRho = 0.95;
gers_opts.verbose   = false;

%% ==================== Main evolutionary loop ====================
if options.verbose, fprintf('Starting evolution...\n'); end

for iGen = 1:ea_params.maxGen
    costs = zeros(ea_params.popSize, 1);
    unstable_count = 0;
    repair_count   = 0;
    
    % --- Fitness Evaluation ---
    for i = 1:ea_params.popSize
        gene = pop{i};
        
        num_links = round(gene(idx_nL));
        act_bin   = logical(gene(idx_act));
        sens_bin  = logical(gene(idx_sens));
        
        % Sparsify the fixed K_dense by strongest elements
        K_sparse = zeros(Nu, Nx);
        keep_idx = K_sorted_idx(1:min(num_links, end));
        K_sparse(keep_idx) = K_dense(keep_idx);
        
        % Apply actuator/sensor masks
        K_sparse(~act_bin, :) = 0;
        K_sparse(:, ~sens_bin) = 0;
        
        % --- 3-phase hard repair for unstable individuals ---
        if options.useGersRepair
            rho_pre = max(abs(eig(A + B * K_sparse)));
            if rho_pre >= 1.0
                [K_sparse, gers_info] = gersgorin_stabilize_K( ...
                    A, B, K_sparse, K_dense, gers_opts);
                if gers_info.success
                    repair_count = repair_count + 1;
                    % Sync genotype to match repaired phenotype
                    gene(idx_nL)   = nnz(K_sparse);
                    gene(idx_act)  = double(any(K_sparse, 2))';
                    gene(idx_sens) = double(any(K_sparse, 1));
                    pop{i} = gene;
                end
            end
        end
        
        % Cost evaluation (returns 1e9 if still unstable after repair)
        J = cost_EA(A, B, Q_fixed, R_fixed, K_sparse, costBM, ...
            ea_params.alpha, max_links, options.gers_weight);
        costs(i) = J;
        
        if J >= 1e9
            unstable_count = unstable_count + 1;
        end
    end
    
    errbuffer(iGen)    = unstable_count;
    repairBuffer(iGen) = repair_count;
    
    % --- Selection: sort by cost ---
    [sortedCosts, sortedIdx] = sort(costs);
    pop   = pop(sortedIdx);
    costs = costs(sortedIdx);
    
    bestCost = sortedCosts(1);
    historyBestCost(iGen) = bestCost;
    historyAvgCost(iGen)  = mean(costs(costs < 1e9));
    
    if options.verbose
        fprintf('Gen %3d: Best=%.4f, Avg=%.4f, Unstable=%d/%d, Repaired=%d\n', ...
            iGen, bestCost, historyAvgCost(iGen), unstable_count, ea_params.popSize, repair_count);
    end
    
    % --- Delta K norm for best individual ---
    bestGene  = pop{1};
    keep_best = K_sorted_idx(1:min(round(bestGene(idx_nL)), end));
    K_sp_best = zeros(Nu, Nx);
    K_sp_best(keep_best) = K_dense(keep_best);
    K_sp_best(~logical(bestGene(idx_act)), :)  = 0;
    K_sp_best(:, ~logical(bestGene(idx_sens))) = 0;
    delta_K_norm(iGen) = norm(K_sp_best - K_dense, 2) / (norm(K_dense, 2) + eps);
    
    % =========================================================
    %  Create next generation: softmax + segment-aware crossover
    % =========================================================
    newPop = cell(ea_params.popSize, 1);
    newPop(1:ea_params.nTop) = pop(1:ea_params.nTop);   % Elitism
    
    % Adaptive softmax temperature
    validCosts = costs(costs < 1e9);
    if numel(validCosts) >= 2
        tau = max(std(validCosts), 1e-6);
    else
        tau = 1.0;
    end
    
    % Softmax selection probabilities (computed once per generation)
    logits = -costs / tau;
    logits(costs >= 1e9) = -inf;
    logits = logits - max(logits);
    weights = exp(logits);
    sumW = sum(weights);
    if sumW < eps
        probs = ones(ea_params.popSize, 1) / ea_params.popSize;
    else
        probs = weights / sumW;
    end
    cumP = cumsum(probs);
    
    for i = (ea_params.nTop + 1):ea_params.popSize
        
        % --- Softmax parent selection ---
        pidx1 = find(cumP >= rand(), 1, 'first');
        pidx2 = find(cumP >= rand(), 1, 'first');
        p1 = pop{pidx1};   J1 = costs(pidx1);
        p2 = pop{pidx2};   J2 = costs(pidx2);
        
        % Fitness advantage weight: w1 -> 1 when parent1 much fitter
        if J1 < 1e9 && J2 < 1e9
            w1 = exp(-J1/tau) / (exp(-J1/tau) + exp(-J2/tau));
        else
            w1 = 0.5;
        end
        
        % --- Segment-aware crossover ---
        child = zeros(1, geneLength);
        
        if rand() >= ea_params.pCross
            if J1 <= J2, child = p1; else, child = p2; end
        else
            % Integer (n_links): fitness-weighted interpolation
            child(idx_nL) = round(w1 * p1(idx_nL) + (1 - w1) * p2(idx_nL));
            
            % Binary (b_act, b_sens): per-bit, P(parent1) = w1
            for j = [idx_act, idx_sens]
                if rand() < w1
                    child(j) = p1(j);
                else
                    child(j) = p2(j);
                end
            end
        end
        
        % --- Mutation ---
        if rand < ea_params.pMutate
            child(idx_nL) = child(idx_nL) + randi([-5, 5]);
        end
        for j = [idx_act, idx_sens]
            if rand < ea_params.pMutate
                child(j) = 1 - child(j);
            end
        end
        
        % --- Bound enforcement ---
        child(idx_nL) = min(max(round(child(idx_nL)), min_links), max_links);
        
        newPop{i} = child;
    end
    
    pop = newPop;
end

%% ==================== Extract results ====================
bestGene       = pop{1};
num_links_best = round(bestGene(idx_nL));
act_bin_best   = logical(bestGene(idx_act));
sens_bin_best  = logical(bestGene(idx_sens));

keep_idx = K_sorted_idx(1:min(num_links_best, end));
K_sparse_best = zeros(Nu, Nx);
K_sparse_best(keep_idx) = K_dense(keep_idx);
K_sparse_best(~act_bin_best, :)  = 0;
K_sparse_best(:, ~sens_bin_best) = 0;

% Final safety net: hard repair on extracted result
if options.useGersRepair
    rho_final = max(abs(eig(A + B * K_sparse_best)));
    if rho_final >= 1.0
        [K_sparse_best, ~] = gersgorin_stabilize_K( ...
            A, B, K_sparse_best, K_dense, gers_opts);
    end
end

% Pack results
result.K_sparse         = K_sparse_best;
result.Q_best           = Q_fixed;
result.R_best           = R_fixed;
result.cost             = historyBestCost(end);
result.gene             = bestGene;
result.active_actuators = find(act_bin_best);
result.active_sensors   = find(sens_bin_best);
result.num_links        = nnz(K_sparse_best);
cost_BM = get_lqr_cost(A, B, Q_fixed, R_fixed, K_bm);
result.costBM       = cost_EA(A, B, Q_fixed, R_fixed, K_bm, cost_BM, ea_params.alpha, max_links);
result.improvement  = (result.costBM - result.cost) / result.costBM;

history.bestCost      = historyBestCost;
history.avgCost       = historyAvgCost;
history.unstableCount = errbuffer;
history.repairCount   = repairBuffer;
history.deltaK_norm   = delta_K_norm;
history.generations   = 1:ea_params.maxGen;

%% ==================== Summary ====================
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
