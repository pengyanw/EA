clear; clc;
addpath(genpath(pwd));  

%% Generate plant
gridSize       = 5; % gridSize by gridSize
connectThresh  = 0.5;
Ts             = 0.2;
actDensity     = 1;
errbuffer = [];
% Want to pick seeds such that grid is fully connected
% (Otherwise decentralized; too easy)
seed = 10; 

% Generate and visualize plants
numNodes    = gridSize*gridSize;
[adjMtx, nodeCoords, susceptMtx, inertiasInv, dampings] = generate_grid_topology(gridSize, connectThresh, seed);

plot_graph(adjMtx, nodeCoords, 'k');

numActs       = round(actDensity*numNodes);
actuatedNodes = randsample(numNodes, numActs);    
sys           = generate_grid_plant(actuatedNodes, adjMtx, susceptMtx, inertiasInv, dampings, Ts);

Nx = sys.Nx;
Nu = sys.Nu;

specRadEA = 1; 
A  = sys.A*specRadEA;
lamda = 0.1;
%A = A + lamda*eye(size(A));
B_ = sys.B2;
%B_= B_ + 1e-1*randn(size(B_));


C_ = eye(Nx);

Q  = eye(Nx);
R_ = eye(sys.Nu);
K_ = -dlqr(A, B_, Q, R_);
suppK_init = zeros(size(K_));
suppK_init(abs(K_)>=1e-2) = true;
K_bar = [];
baseCost = get_lqr_cost(A, B_, Q, R_, K_);
alpha = 0.9; %select those K s.t. eigvals larger than 5 times basecost
lqr_thresh = alpha*baseCost;
rank(B_)==size(B_,2)
checkctrb(A,B_)
% ===== check if (A,B) (L, alpha)-stable =====

spectral_radius = max(abs(eig(A)));
is_alpha_stable = spectral_radius <= alpha;
new = A+B_*K_;
new_spectral_radius = max(abs(eig(new)))




%% EA
% Hyperparameters
popSize = 20;     % Population size
pMutate = 0.1;
nTop    = 15;      % Historical "best" individuals
genDiffStop = 1e2; % If we go this many generations without improvement then stop

sensPen = 0.2; % one sensor costs same as Z% of LQR degradation
actPen  = 0.4; % one actuator costs same as X% of LQR degradation
commPen = 0.05; % one communication link costs same as Y% of LQR degradation 

currPopBin = cell(popSize, 1); % which comm links to use
currPopCts = cell(popSize, 1); % [# sensors, # actuators]
bestGens   = zeros(nTop, 1);
bestCosts  = inf(nTop, 1);
bestChromsBin = cell(nTop, 1);
bestChromsCts = cell(nTop, 1);

% Initialization
currGen = 1;
%record every gen's best individual(lowest cost)
historyGen = [];
historyBestCost = [];
costlqr = get_lqr_cost(A,B_,Q,R_,K_);    
for i = 1:popSize
    if i <= min(ceil(popSize/length(popSize))-1,1)
        % First half: same sparsity pattern as K_
        currPopBin{i} = reshape(logical(suppK_init), 1, []);
        currPopCts{i} = [(Nx), (Nu)];
   else
        % Last half: Random Bin
        bin = logical(randi(2, 1, Nx * Nu));
        currPopBin{i} = bin;

        % === [New] Use row/column L1 scores to pick important actuators/sensors ===
        K_mask = reshape(bin, Nu, Nx);
        row_strength = sum(abs(K_ .* K_mask), 2);  % actuator importance
        col_strength = sum(abs(K_ .* K_mask), 1);  % sensor importance

        % Select number of sensors/actuators based on percentiles or random top-k
        nKeepRow = randi(Nu);  % e.g., random number to keep
        nKeepCol = randi(Nx);

        [~, row_idx] = sort(row_strength, 'descend');
        [~, col_idx] = sort(col_strength, 'descend');

        activeRows = row_idx(1:nKeepRow);
        activeCols = col_idx(1:nKeepCol);

        % Final counts of retained actuator/sensor channels
        currPopCts{i} = [length(activeCols), length(activeRows)];
    end
end
    
    


lastImproved = currGen; % When were the top individuals last improved?
%% evolution

tic;
while true
    if (currGen - lastImproved) >= genDiffStop
        fprintf('\nNo changes to the top %d have happened in %d generations\n', ...
                    nTop, genDiffStop);
        fprintf('Stopping at generation %d.\n', currGen);
        break
    end
    %[updatedBin, Kstar] = sparsify_lqr_fsolve(A, B_, Q, R_, K_, currPopBin, alpha);
    [costs, currK] = cost(currPopBin, currPopCts, A, B_, C_, Q, R_, K_, baseCost, sensPen, actPen, commPen);

    if ~exist('unstable','var'); unstable = {}; end
    unstable{currGen} = struct('bin', {}, 'cts', {}, 'cost', {}, 'K', {});

    for i = 1:length(currPopBin)
        if ~isfinite(costs(i)) || costs(i) >= 1e4*baseCost
            unstable{currGen}(end+1) = struct( ...
            'bin',  currPopBin{i}, ...   % 1 x (Nx*Nu) logical row vector
            'cts',  currPopCts{i}, ...   % [#sensors, #actuators]
            'cost', costs(i),...           % store the cost as well
            'K', currK(i,:,:));
        end
    end
%     [costs, currK] = cost_nocts(currPopBin, currPopCts, A, B_, C_, Q, R_, K_, baseCost, sensPen, actPen, commPen);
    K_bar(currGen,:,:,:) = currK;

    if max(costs) == inf
    errbuffer(currGen) = length(find(costs==inf));
    costs(find(costs==inf))=max(costs(costs~=inf))*1e6;
    
    end

    for i=1:popSize % Fitness/cost evaluation
        idx      = find(bestCosts > costs(i), 1, 'first');%select highest cost population 
        if ~isempty(idx) % Add to top list
            lastImproved = currGen;
            bestGens   = [bestGens(1:idx-1); currGen; bestGens(idx:end-1)];
            bestCosts  = [bestCosts(1:idx-1); costs(i); bestCosts(idx:end-1)];
            bestChromsBin = {bestChromsBin{1:idx-1} currPopBin{i} bestChromsBin{idx:end-1}};
            bestChromsCts = {bestChromsCts{1:idx-1} currPopCts{i} bestChromsCts{idx:end-1}};
        end
    end
    
    if lastImproved == currGen % Only print if you've updated something
        print_top_list(currGen, bestGens, bestCosts, bestChromsCts)
%         historybestBin = currPopBin;
    end
    
    repProb = generate_reproduction_probs(costs,costlqr);
    [currPopBin, currPopCts] = generate_next_pop(currPopBin, currPopCts, repProb, pMutate, Nx, Nu);
    currGen = currGen+1;
    historyGen(end+1) = currGen;
    historyBestCost(end+1) = bestCosts(1); % current lowest cost individual

end
toc



%% Evaluation 1
% Benchmark: semi-truncated K
KSuppBM     = abs(K_) > 1e-2;
K1          = zeros(size(K_));
K1(KSuppBM) = K_(KSuppBM);

specRadBM  = max(abs(eig(A+B_*K1)));
lqrCostBM  = get_lqr_cost(A, B_, Q, R_, K1) ./ baseCost;
costBM     = lqrCostBM + sensPen*Nx + actPen*Nu + commPen*nnz(K1);

fprintf('\n==========Benchmark:==========\n')
fprintf('stab: %.2f, lqr: %f, total: %f\n', specRadBM, lqrCostBM, costBM)
fprintf('sens: %d, act: %d, comm: %d\n', Nx, Nu, nnz(K1));

% Best option from EA
[YIdx, UIdx] = postprocess(bestChromsBin, bestChromsCts, K_);

bin = bestChromsBin{1}; cts = bestChromsCts{1};

KSupp      = reshape(bin, size(K_));
K2         = K_;
K2(~KSupp) = 0;
K2(:, ~YIdx{1}) = []; % Eliminate unsensed columns
K2(~UIdx{1}, :) = []; % Eliminates unactuated rows

B = B_(:, UIdx{1});
C = C_(YIdx{1}, :);
R = R_(UIdx{1}, UIdx{1});

specRadEA  = max(abs(eig(A+B*K2*C)));
lqrCostEA  = get_lqr_cost(A, B, Q, R, K2*C) ./ baseCost;
costEA     = lqrCostEA + sensPen*cts(1) + actPen*cts(2) + commPen*nnz(K2);

fprintf('\n==========EA output:==========\n')
fprintf('stab: %.2f, lqr: %f, total: %f\n', specRadEA, lqrCostEA, costEA)
fprintf('sens: %d, act: %d, comm: %d\n', cts(1), cts(2), nnz(K2));

fprintf('\nImprovement over benchmark: %.2f\n', (costBM - costEA)/costBM)

%% New LQR based on EA
Knew = -dlqr(A, B, Q, R);
K3   = zeros(size(K_));
K3(UIdx{1}, :) = Knew;
K3(~KSupp) = 0;
K3(:, ~YIdx{1}) = []; % Eliminate unsensed columns
K3(~UIdx{1}, :) = []; % Eliminates unactuated rows

specRadEA2  = max(abs(eig(A+B*K3*C)));
lqrCostEA2  = get_lqr_cost(A, B, Q, R, K3*C) ./ baseCost;
costEA2     = lqrCostEA2 + sensPen*cts(1) + actPen*cts(2) + commPen*nnz(K3);

fprintf('\n==========EA-based LQR:==========\n')
fprintf('stab: %.2f, lqr: %f, total: %f\n', specRadEA2, lqrCostEA2, costEA2)
fprintf('sens: %d, act: %d, comm: %d\n', cts(1), cts(2), nnz(K3));
disp(K3)
fprintf('\nImprovement over benchmark: %.2f\n', (costBM - costEA2)/costBM)
%% 
% 
figure;
plot(historyGen, historyBestCost, 'r-o', 'LineWidth', 1.5);
xlabel('Generation');
ylabel('Best Cost');
title('Evolution of Best Cost per Generation');
grid on;

% 
[minCost, minIdx] = min(historyBestCost);
fprintf('\n==========EA Summary==========\n');
fprintf('Best cost: %f at generation %d\n', minCost, historyGen(minIdx));
% save figure
if ~exist('figures_updatedCts', 'dir')
    mkdir('figures_updatedCts');
end

improvementRatio = (costBM - costEA) / costBM;

fileName = sprintf(['figures_updatedCts/evolution_g%d_p%.2f_pop%d_top%d_stop%d_' ...
                    'alpha%.2f_advVSbechmark%.2f.png'], ...
                    gridSize, pMutate, popSize, nTop, genDiffStop, ...
                    alpha, improvementRatio);



saveas(gcf, fileName);
fprintf('Saved evolution plot to %s\n', fileName);
%% Plot members that have infinite LQR cost in every generation
%  errbuffer: 1 x N double, representing times Inf appear

gens = 1:length(errbuffer);
figure;
plot(gens, errbuffer, 'o-', 'LineWidth', 2);
xlabel('Generation');
ylabel('# Inf LQR Costs');
title(sprintf('Inf LQR Cost per Generation (grid=%d, pop=%d, top=%d, inf_cnt=%d)', gridSize, popSize, nTop, mean(errbuffer)));
grid on;
ylim([0 popSize]);
fprintf('mean unstable cnts: %.2f\n', mean(errbuffer));

% Save
folder = 'Inf buffer';
if ~exist(folder, 'dir'); mkdir(folder); end
filename = sprintf('%s/infLQR_grid%d_pop%d_top%d.png', folder, gridSize, popSize, nTop);
saveas(gcf, filename);
fprintf('Figure saved to: %s\n', filename);
save(sprintf('unstable_grid%d_pop%d_top%d.mat', gridSize, popSize, nTop), ...
     'unstable', 'errbuffer', 'gridSize', 'popSize', 'nTop');

