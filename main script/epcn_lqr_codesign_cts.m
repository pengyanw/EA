clear; clc;
addpath(genpath(pwd));  
close all
%% 1. System Generation (No changes here)
% =========================================================================
gridSize       = 5;
connectThresh  = 0.5;
Ts             = 0.2;
actDensity     = 1;
seed           = 9; 

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
%% 2. EA Parameters (Updated for the new strategy)
% =========================================================================
popSize     = 20;      % Population size
maxGen      = 150;      % Maximum generations
pMutate     = 0.05;      % Mutation probability for each gene component
pCross      = 0.8;      % Crossover probability
nTop        = 10;       % Number of top individuals (elites) to keep
alpha       = 0;      % Cost function weighting factor (LQR vs hardware)

% --- NEW: Parameters for the new gene structure ---
% The gene is now: [diag(Q) elements, diag(R) elements, num_links]
% --- NEW: Extended Gene Structure ---
len_diag_Q    = Nx;
len_diag_R    = Nu;
len_num_links = 1;
len_act_bin   = Nu;
len_sens_bin  = Nx;
geneLength    = len_diag_Q + len_diag_R + len_num_links + len_act_bin + len_sens_bin;


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
    % Part 1: Diagonal Q
    diag_Q = (max_Q_val - min_Q_val) * rand(1, len_diag_Q) + min_Q_val;
    
    % Part 2: Diagonal R
    diag_R = (max_R_val - min_R_val) * rand(1, len_diag_R) + min_R_val;
    
    % Part 3: Communication link count
    num_links = randi([min_links, max_links]);
    
    % Part 4: Actuator and Sensor binary masks
    act_bin  = double(rand(1, len_act_bin)  > 0.01);
    sens_bin = double(rand(1, len_sens_bin) > 0.01);

    % Combine into full gene
    pop{i} = [diag_Q, diag_R, num_links, act_bin, sens_bin];
end


% History tracking
historyBestCost = zeros(maxGen, 1);
historyAvgCost  = zeros(maxGen, 1);
historyGen      = 1:maxGen;
errbuffer       = zeros(maxGen, 1); % Tracks number of unstable individuals
popbuffer = [];
flag=  []; %validate the theoretical bound 
%% 5. Main Evolutionary Algorithm Loop (Updated with bin+cts structure)
% =========================================================================
fprintf('Starting evolution...\n');
delta_K_norm = []; % Record K change ratio

for iGen = 1:maxGen
    costs = zeros(popSize, 1);
    unstable_count = 0;
    cts = cell(popSize, 1);  % store [active_rows, active_cols] for each individual
    
    % --- Fitness Evaluation ---
    for i = 1:popSize
        gene = pop{i};

        % Step A: Decode gene parts
        diag_Q    = gene(1:len_diag_Q);
        diag_R    = gene(len_diag_Q + 1 : len_diag_Q + len_diag_R);
        num_links = round(gene(len_diag_Q + len_diag_R + 1));

        act_bin   = logical(gene(len_diag_Q + len_diag_R + len_num_links + (1:len_act_bin)));
        sens_bin  = logical(gene(end - len_sens_bin + 1:end));

        Q = diag(diag_Q);
        R = diag(diag_R);

        % Step B: Compute dense controller
        try
            K_dense = -dlqr(A, B_, Q, R);
        catch
            costs(i) = 1e12;
            unstable_count = unstable_count + 1;
            cts{i} = [0 0];
            continue;
        end

        % Step C: Sparsify K_dense by strongest |K| elements
        K_sparse = zeros(Nu, Nx);
        [~, sorted_idx] = sort(abs(K_dense(:)), 'descend');
        keep_idx = sorted_idx(1:min(num_links, end));
        K_sparse(keep_idx) = K_dense(keep_idx);

        % Step D: Apply actuator/sensor binary masks
        K_sparse(~act_bin, :) = 0;   % deactivate rows (actuators)
        K_sparse(:, ~sens_bin) = 0;  % deactivate cols (sensors)

        % Step E: Compute cost using external cost_EA()
        J = cost_EA(A, B_, Q_bm, R_bm_, K_sparse, costBM, alpha, max_links);
        costs(i) = J;
        flag = [flag check_stability_margin(A,B_,K_sparse,K_dense)];
        % Step F: Count active structure
        active_rows = nnz(any(K_sparse,2));
        active_cols = nnz(any(K_sparse,1));
        cts{i} = [active_rows, active_cols];

        if J >= 1e9
            unstable_count = unstable_count + 1;
        end
    end

    errbuffer(iGen) = unstable_count;
    %save(fullfile('cache', sprintf('cts_gen_%03d.mat', iGen)), 'cts');

    % --- Selection, Crossover, and Mutation ---
    [sortedCosts, sortedIdx] = sort(costs);
    pop = pop(sortedIdx); % sort by fitness
    popbuffer = [popbuffer pop];

    bestCost = sortedCosts(1);
    historyBestCost(iGen) = bestCost;
    historyAvgCost(iGen)  = mean(costs(costs < 1e9));
    fprintf('Gen %d: Best Cost=%.4f, Avg=%.4f, Unstable=%d/%d\n', ...
        iGen, bestCost, historyAvgCost(iGen), unstable_count, popSize);

    % Compute best controller of this generation
    bestGene = pop{1};
    diag_Q_best = bestGene(1:len_diag_Q);
    diag_R_best = bestGene(len_diag_Q + 1 : len_diag_Q + len_diag_R);
    Q_best = diag(diag_Q_best);
    R_best = diag(diag_R_best);
    best_K_dense = -dlqr(A, B_, Q_best, R_best);
    [~, sorted_idx] = sort(abs(best_K_dense(:)), 'descend');
    keep_idx = sorted_idx(1:min(bestGene(len_diag_Q + len_diag_R + 1), end));
    K_sparse_best = zeros(Nu, Nx);
    K_sparse_best(keep_idx) = best_K_dense(keep_idx);
    K_sparse_best(~logical(bestGene(len_diag_Q + len_diag_R + len_num_links + (1:len_act_bin))), :) = 0;
    K_sparse_best(:, ~logical(bestGene(end - len_sens_bin + 1:end))) = 0;

    delta_K_norm(iGen) = norm(K_sparse_best-best_K_dense, 2) / (norm(best_K_dense, 2) + eps);

    % --- Create the next generation ---
    newPop = cell(popSize, 1);
    newPop(1:nTop) = pop(1:nTop); % elitism

    for i = (nTop+1):popSize
        % Parent selection
        parent1 = pop{randi(popSize)};
        parent2 = pop{randi(popSize)};

        % Crossover
        if rand < pCross
            crossPoint = randi(geneLength - 1);
            child = [parent1(1:crossPoint), parent2(crossPoint+1:end)];
        else
            child = parent1;
        end

        % Mutation (continuous parts)
        for j = 1:(len_diag_Q + len_diag_R + len_num_links)
            if rand < pMutate
                child(j) = child(j) + randn * 0.1 * (max_Q_val - min_Q_val);
            end
        end
        
        % Bound enforcement
        child(1:len_diag_Q) = min(max(child(1:len_diag_Q), min_Q_val), max_Q_val);
        child(len_diag_Q+1 : len_diag_Q+len_diag_R) = ...
            min(max(child(len_diag_Q+1 : len_diag_Q+len_diag_R), min_R_val), max_R_val);
        child(len_diag_Q+len_diag_R+1) = ...
            min(max(round(child(len_diag_Q+len_diag_R+1)), min_links), max_links);

        % Mutation for binary parts
        for j = (len_diag_Q + len_diag_R + len_num_links + 1):geneLength
            if rand < pMutate
                child(j) = 1 - child(j);
            end
        end

        newPop{i} = child;
    end

    pop = newPop;
end



%% 6. Results and Plotting (Updated for clarity)
close all
% =========================================================================
[minCost, minIdx] = min(historyBestCost);
bestGen = historyGen(minIdx);
costEA = minCost; % Final cost from EA

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
saveas(figure(1), sprintf('figures/evo_bestcost_grid%dseed%d.png', gridSize, seed));
saveas(figure(2), sprintf('figures/evo_avgcost_grid%dseed%d.png', gridSize, seed));
saveas(figure(3), sprintf('figures/unstable_count_grid%dseed%d.png', gridSize, seed));
fprintf("max deltaK/K condition number=%d", max(delta_K_norm))
[res,~,~] = trim_nonzero(K_sparse_best);


save(fullfile('cache', 'popbuffer.mat'), 'popbuffer');