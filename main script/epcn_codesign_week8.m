%% Evolutionary Algorithm with Binary Mask Gene
clear; clc;
addpath(genpath(pwd));  
close all

%% 1. System Generation
% =========================================================================
gridSize       = 5;
connectThresh  = 0.5;
Ts             = 0.2;
actDensity     = 1;
seed           = 5; 

numNodes    = gridSize*gridSize;
[adjMtx, nodeCoords, susceptMtx, inertiasInv, dampings] = ...
    generate_grid_topology(gridSize, connectThresh, seed);
plot_graph(adjMtx, nodeCoords, 'k');

numActs       = round(actDensity*numNodes);
actuatedNodes = randsample(numNodes, numActs);    
sys           = generate_grid_plant(actuatedNodes, adjMtx, ...
    susceptMtx, inertiasInv, dampings, Ts);

Nx = sys.Nx;
Nu = sys.Nu;
A  = sys.A;
B_ = sys.B2;

%% 2. EA Parameters
% =========================================================================
popSize     = 20;      % Population size
maxGen      = 150;     % Maximum generations
pMutate     = 0.1;     % Mutation probability
pCross      = 0.8;     % Crossover probability
nTop        = 10;      % Elitism (number of elites)
alpha       = 0.5;     % Weighting between LQR and sparsity cost

% --- Gene Structure ---
len_diag_Q = Nx;
len_diag_R = Nu;
len_mask   = Nu * Nx;   % Binary mask for sparsity
geneLength = len_diag_Q + len_diag_R + len_mask;

% --- Value Bounds ---
max_Q_val = 10; min_Q_val = 1e-4;
max_R_val = 10; min_R_val = 1e-4;

%% 3. Benchmark Calculation (No changes here)
% =========================================================================
Q_bm  = eye(Nx);
R_bm_ = eye(Nu);
K_bm  = -dlqr(A, B_, Q_bm, R_bm_); % Benchmark dense LQR controller
costBM = get_lqr_cost(A, B_, Q_bm, R_bm_, K_bm);
fprintf('Benchmark LQR cost (dense controller): %f\n', costBM);

% Benchmark: semi-truncated K
KSuppBM     = abs(K_bm) > 1e-2;
K1          = zeros(size(K_bm));
K1(KSuppBM) = K_bm(KSuppBM);
cost_bm     = cost_EA(A, B_, Q_bm, R_bm_, K1, costBM, alpha)

%% 4. Initialization
% =========================================================================
fprintf('Initializing population (binary mask)...\n');
pop = cell(popSize, 1);
for i = 1:popSize
    diag_Q = (max_Q_val - min_Q_val) * rand(1, len_diag_Q) + min_Q_val;
    diag_R = (max_R_val - min_R_val) * rand(1, len_diag_R) + min_R_val;
    mask_vec = rand(1, len_mask) < 0.2; % Initial sparsity rate (20%)
    pop{i} = [diag_Q diag_R 1.*(reshape(K1, 1, [])>0)];   % 转为 1 × (Nu*Nx) 行向量
end

% History variables
historyBestCost = zeros(maxGen, 1);
historyAvgCost  = zeros(maxGen, 1);
errbuffer       = zeros(maxGen, 1);
popbuffer       = {};
flag = [];

%% 5. Main Evolutionary Algorithm Loop
% =========================================================================
fprintf('Starting evolution...\n');
for iGen = 1:maxGen
    costs = zeros(popSize, 1);
    unstable_count = 0;
    
    for i = 1:popSize
        gene = pop{i};
        diag_Q = gene(1:len_diag_Q);
        diag_R = gene(len_diag_Q + 1 : len_diag_Q + len_diag_R);
        mask_vec = gene(len_diag_Q + len_diag_R + 1 : end);
        mask = reshape(mask_vec, [Nu, Nx]);
        
        Q = diag(diag_Q);
        R = diag(diag_R);
        
        % Compute dense LQR controller
        try
            K_dense = -dlqr(A, B_, Q, R);
        catch
            costs(i) = 1e12; 
            unstable_count = unstable_count + 1;
            continue;
        end
        
        % Apply evolved mask
        K_sparse = K_dense .* mask;
        
        % Check stability
        if max(abs(eig(A + B_ * K_sparse))) >= 1.0
            costs(i) = 1e9;
            unstable_count = unstable_count + 1;
        else
            costs(i) = cost_EA(A, B_, Q_bm, R_bm_, ...
                K_sparse, costBM, alpha, Nu*Nx);
        
        end
        flag = [flag check_stability_margin(A,B_,K_sparse,K_dense)];
    end
    
    errbuffer(iGen) = unstable_count;
    
    % --- Selection ---
    [sortedCosts, sortedIdx] = sort(costs);
    pop = pop(sortedIdx);
    popbuffer = [popbuffer pop];
    bestCost = sortedCosts(1);
    
    historyBestCost(iGen) = bestCost;
    historyAvgCost(iGen)  = mean(costs(costs < 1e9));
    
    fprintf('Gen %d: Best=%.4f | Avg=%.4f | Unstable=%d/%d\n', ...
        iGen, bestCost, historyAvgCost(iGen), unstable_count, popSize);
    
    % --- Reproduction ---
    newPop = cell(popSize,1);
    newPop(1:nTop) = pop(1:nTop);
    
    for i = (nTop+1):popSize
        % Parent selection
        parent1 = pop{randi(popSize)};
        parent2 = pop{randi(popSize)};
        
        % Crossover
        if rand < pCross
            crossPoint = randi(geneLength-1);
            child = [parent1(1:crossPoint), parent2(crossPoint+1:end)];
        else
            child = parent1;
        end
        
        % Mutation (Q and R)
        for j = 1:(len_diag_Q + len_diag_R)
            if rand < pMutate
                child(j) = child(j) + randn * 0.1 * (max_Q_val - min_Q_val);
                child(j) = max(min(child(j), max_Q_val), min_Q_val);
            end
        end
        
        % Mutation (binary mask flipping)
        mask_start = len_diag_Q + len_diag_R + 1;
        mask_end   = geneLength;
        for j = mask_start:mask_end
            if rand < pMutate
                child(j) = 1 - child(j);
            end
        end
        
        newPop{i} = child;
    end
    pop = newPop;
end

%% 6. Results and Plotting
% =========================================================================
[minCost, minIdx] = min(historyBestCost);
bestGen = minIdx;
bestGene = pop{1};
mask_best = reshape(bestGene(len_diag_Q+len_diag_R+1:end), [Nu, Nx]);
Q_best = diag(bestGene(1:len_diag_Q));
R_best = diag(bestGene(len_diag_Q+1 : len_diag_Q+len_diag_R));
K_dense_best = -dlqr(A,B_,Q_best,R_best);
K_sparse_best = K_dense_best .* mask_best;

fprintf('\n========== EA Summary ==========\n');
fprintf('Best cost found: %.4f at generation %d\n', minCost, bestGen);
fprintf('Nonzero links: %d / %d\n', nnz(mask_best), numel(mask_best));

figure;
subplot(1,2,1);
imagesc(mask_best);
title('Evolved sparsity mask');
xlabel('State index'); ylabel('Actuator index');
colormap(gray);

subplot(1,2,2);
plot(1:maxGen, historyBestCost, 'b-', 'LineWidth', 2); hold on;
plot(1:maxGen, historyAvgCost, 'r--', 'LineWidth', 1.5);
xlabel('Generation'); ylabel('Cost');
legend('Best', 'Average');
title('EA Cost Evolution'); grid on;

figure;
plot(1:maxGen, errbuffer, 'r-o', 'LineWidth', 1.5);
xlabel('Generation'); ylabel('Unstable Individuals');
title('Instability Count per Generation'); grid on;

% Save
if ~exist('cache','dir'), mkdir('cache'); end
save(fullfile('cache','EA_binary_mask_results.mat'), 'popbuffer', 'mask_best', 'K_sparse_best');
