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
A = A + lamda*eye(size(A));
B_ = sys.B2;
%B_= B_ + 1e-1*randn(size(B_));
%% compute the graph properties
% min distance
G = graph(adjMtx);
D = distances(G);           % Nnode x Nnode, D(i,j)=最短跳数
D(isinf(D)) = max(D(~isinf(D)))+1;  % 断开分量用大值兜底
n_s = 2;                          % 每格点状态数（按你的模型）

n_s  = Nx/ Nu;
D_for_K = kron(D, ones(1, n_s));  % Nu x Nx
use_gaussian = false;   % false 则用 exp(-beta*d)
sigma  = 1.5;          % 高斯核带宽
beta   = 0.7;          % 指数核衰减系数（备选）
w_min  = 0.0;          % 门控下限（可设 0 ~ 0.2 之间）

if use_gaussian
    W = exp(-(D_for_K./sigma).^2);
else
    W = exp(-beta * D_for_K);
end
W = max(W, w_min);


%% 2. EA Parameters (Updated for the new strategy)
% =========================================================================
popSize     = 20;      % Population size
maxGen      = 150;      % Maximum generations
pMutate     = 0.1;      % Mutation probability for each gene component
pCross      = 0.8;      % Crossover probability
nTop        = 10;       % Number of top individuals (elites) to keep
alpha       = 0;      % Cost function weighting factor (LQR vs hardware)

% --- NEW: Parameters for the new gene structure ---
% The gene is now: [diag(Q) elements, diag(R) elements, num_links]
len_diag_Q    = Nx;
len_diag_R    = Nu;
len_num_links = 1;
geneLength    = len_diag_Q + len_diag_R + len_num_links;

% --- NEW: Bounds for gene values ---
max_Q_val   = 10;   % Max value for a diagonal element of Q
max_R_val   = 10;   % Max value for a diagonal element of R
min_Q_val   = 1e-4; % Min value to keep matrices positive semi-definite
min_R_val   = 1e-4; % Min value to keep matrices positive definite
min_links   = 1;    % Minimum number of non-zero elements in K
max_links   = Nu * Nx; % Maximum number of non-zero elements in K

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
%% 4. EA Initialization (Completely Rewritten)
% =========================================================================
fprintf('Initializing population...\n');
pop = cell(popSize, 1);
for i = 1:popSize
    % Gene part 1: Diagonal elements of Q
    diag_Q = (max_Q_val - min_Q_val) * rand(1, len_diag_Q) + min_Q_val;
    
    % Gene part 2: Diagonal elements of R
    diag_R = (max_R_val - min_R_val) * rand(1, len_diag_R) + min_R_val;
    
    % Gene part 3: Number of communication links
    num_links = randi([min_links, nnz(KSuppBM)]);
    
    % Combine to form the chromosome
    pop{i} = [diag_Q, diag_R, num_links];
end

% History tracking
historyBestCost = zeros(maxGen, 1);
historyAvgCost  = zeros(maxGen, 1);
historyGen      = 1:maxGen;
errbuffer       = zeros(maxGen, 1); % Tracks number of unstable individuals
popbuffer = [];
flag = [];
%% 5. Main Evolutionary Algorithm Loop (Heavily Modified)
% =========================================================================
fprintf('Starting evolution...\n');
delta_K_norm = []; %record evry round's K norm/ K_
for iGen = 1:maxGen
    
    costs = zeros(popSize, 1);
    unstable_count = 0;
    
    % --- Fitness Evaluation ---
    for i = 1:popSize
        gene = pop{i};
        
        % Step A: Decode gene into Q, R, and num_links
        diag_Q    = gene(1:len_diag_Q);
        diag_R    = gene(len_diag_Q + 1 : len_diag_Q + len_diag_R);
        num_links = round(gene(end)); % Ensure it's an integer
        
        Q = diag(diag_Q);
        R = diag(diag_R);
        
        % Step B: Solve for the dense, stable controller K_dense
        % Use a try-catch block in case dlqr fails (e.g., if R is not pos-def)
        try
            K_dense = -dlqr(A, B_, Q, R);
        catch
            costs(i) = 1e12; % Assign a massive penalty if dlqr fails
            unstable_count = unstable_count + 1;
            continue;
        end
        
        % Step C: Sparsify K_dense
        K_sparse = zeros(Nu, Nx);
        [~, sorted_idx] = sort(abs(K_dense(:)), 'descend');
        keep_indices = sorted_idx(1:min(num_links, end));
        K_sparse(keep_indices) = K_dense(keep_indices);
        
        % Step D: Check stability of the sparse controller
        if max(abs(eig(A + B_ * K_sparse))) >= 1.0
            costs(i) = 1e9; % Assign a large penalty for instability
            unstable_count = unstable_count + 1;
        else
            
            % Normalized and weighted final cost
            costs(i) = cost_EA(A, B_, Q_bm, R_bm_, K_sparse, costBM, alpha, max_links);
            flag = [flag check_stability_margin(A,B_,K_sparse,K_dense)];
        end
    end
    
    errbuffer(iGen) = unstable_count;
    
    % --- Selection, Crossover, and Mutation ---
    [sortedCosts, sortedIdx] = sort(costs);
    pop = pop(sortedIdx); % Sort population by fitness
    popbuffer = [popbuffer pop];
    best_K = -dlqr(A, B_, diag(pop{1}(1:len_diag_Q)), diag(pop{1}(len_diag_Q + 1 : len_diag_Q + len_diag_R)));
    [~, sorted_idx] = sort(abs(best_K(:)), 'descend');
    keep_indices = sorted_idx(1:min(pop{1}(end), end));
    K_sparse_best = zeros(Nu, Nx);
    K_sparse_best(keep_indices) = best_K(keep_indices);
    
    delta_K_norm(iGen) = norm(K_sparse-best_K,2)/norm(best_K,2);
    bestIndividual = [pop{1} ];
    bestCost = sortedCosts(1);
    
    historyBestCost(iGen) = bestCost;
    historyAvgCost(iGen)  = mean(costs(costs < 1e9)); % Avg cost of stable solutions
    
    fprintf('Gen %d: Best Cost=%.4f, Avg Cost=%.4f, Unstable=%d/%d\n', ...
            iGen, bestCost, historyAvgCost(iGen), unstable_count, popSize);
    
    % Create the next generation
    newPop = cell(popSize, 1);
    
    % Elitism: Keep the top individuals
    newPop(1:nTop) = pop(1:nTop);
    
    % Fill the rest of the population
    for i = (nTop + 1):popSize
        % Select parents (e.g., tournament selection or roulette wheel)
        parent1 = pop{randi(popSize)}; 
        parent2 = pop{randi(popSize)};
        
        % Crossover
        if rand < pCross
            crossPoint = randi(geneLength - 1);
            child = [parent1(1:crossPoint), parent2(crossPoint+1:end)];
        else
            child = parent1;
        end
        
        % Mutation
        for j = 1:(geneLength - 1) % Mutate Q and R elements
            if rand < pMutate
                % Add Gaussian noise
                child(j) = child(j) + randn * 0.1 * (max_Q_val - min_Q_val);
                % Enforce bounds
                child(j) = max(min(child(j), max_Q_val), min_Q_val);
            end
        end
        % Mutate num_links separately
        if rand < pMutate
            child(end) = child(end) + randi([-5, 5]); % Small integer change
            % Enforce bounds
            child(end) = max(min(child(end), max_links), min_links);
        end
        
        newPop{i} = child;
    end
    
    pop = newPop;
end

%% 6. Results and Plotting (Updated for clarity)
% =========================================================================
[minCost, minIdx] = min(historyBestCost);
bestGen = historyGen(minIdx);
costEA = minCost; % Final cost from EA
[res,~,~] = trim_nonzero(K_sparse_best);
fprintf("Final K:%d\n", res)
fprintf('\n========== EA Summary ==========\n');
fprintf('Best cost found: %f at generation %d\n', costEA, bestGen);
fprintf('best adavntage over benchmark:%d\n', (-costEA+cost_bm)/cost_bm)
figure;
plot(historyGen, historyBestCost, 'b-', 'LineWidth', 2, 'DisplayName', 'Best Cost');

hold on;

legend;
xlabel('Generation');
ylabel('Normalized Cost');

title('Evolution of LQR Co-Design Cost');
grid on;

figure;
plot(historyGen, historyAvgCost, 'g--', 'LineWidth', 1.5, 'DisplayName', 'Average Stable Cost');
legend;
xlabel('Generation');
ylabel('Normalized Cost');

title('Evolution of average LQR Co-Design Cost');
grid on;

% Plot number of unstable individuals per generation
figure;
plot(historyGen, errbuffer, 'r-o', 'LineWidth', 2);
xlabel('Generation');
ylabel('Number of Unstable Individuals');
title('Stability Failures per Generation');
grid on;
ylim([0 popSize]);

% Ensure output folder exists
if ~exist('figures', 'dir'), mkdir('figures'); end

% Save all current figures with gridSize in filename
saveas(figure(2), sprintf('figures/evo_bestcost_grid%dseed%d.png', gridSize, seed));
saveas(figure(3), sprintf('figures/evo_avgcost_grid%dseed%d.png', gridSize, seed));
saveas(figure(4), sprintf('figures/unstable_count_grid%dseed%d.png', gridSize, seed));
fprintf("max deltaK/K condition number=%d", max(delta_K_norm))


save(fullfile('cache', 'popbuffer.mat'), 'popbuffer');