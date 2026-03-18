%% Multi-Seed A/B Experiment: EA with vs without Gershgorin Repair
% Compare convergence and stability across seeds, report seed-averaged results.
%
% Author: Pengyang Wu
% Date: 2026-02-03

clear; clc; close all;
addpath(genpath(pwd));

%% ===================== Configuration =====================
gridSize       = 4;
connectThresh  = 0.5;
Ts             = 0.2;
actDensity     = 1;
topoSeed       = 10;

eaSeeds  = [1, 5, 10, 17, 42, 66, 88, 100];
numSeeds = length(eaSeeds);

ea_params.popSize  = 20;
ea_params.maxGen   = 120;
ea_params.pMutate  = 0.05;
ea_params.pCross   = 0.8;
ea_params.nTop     = 10;
ea_params.alpha    = 0;

opts_base.max_Q_val  = 10;
opts_base.max_R_val  = 10;
opts_base.min_Q_val  = 1e-4;
opts_base.min_R_val  = 1e-4;
opts_base.verbose    = false;

convergeTol = 0.05;
maxGen = ea_params.maxGen;

%% ===================== System Generation =====================
fprintf('Generating system (grid %dx%d, topoSeed=%d)...\n', gridSize, gridSize, topoSeed);
numNodes = gridSize * gridSize;
[adjMtx, nodeCoords, susceptMtx, inertiasInv, dampings] = ...
    generate_grid_topology(gridSize, connectThresh, topoSeed);
rng(topoSeed);
numActs      = round(actDensity * numNodes);
actuatedNodes = randsample(numNodes, numActs);
sys = generate_grid_plant(actuatedNodes, adjMtx, susceptMtx, inertiasInv, dampings, Ts);
A = sys.A;  B = sys.B2;
[Nx, Nu] = size(B);
fprintf('System: Nx=%d, Nu=%d\n\n', Nx, Nu);

%% ===================== Pre-allocate =====================
% Two conditions: (1) no repair, (2) with Gershgorin repair
condNames = {'Baseline (no repair)', 'Gershgorin repair'};
nCond = 2;

bestCostAll  = zeros(maxGen, numSeeds, nCond);
unstableAll  = zeros(maxGen, numSeeds, nCond);
repairAll    = zeros(maxGen, numSeeds, nCond);
finalCostAll = zeros(numSeeds, nCond);
convGenAll   = zeros(numSeeds, nCond);
convRateAll  = zeros(numSeeds, nCond);
elapsedAll   = zeros(numSeeds, nCond);

%% ===================== Run Experiments =====================
for c = 1:nCond
    opts = opts_base;
    opts.useGersRepair = (c == 2);
    tag = condNames{c};
    fprintf('========== Condition %d: %s ==========\n', c, tag);

    for s = 1:numSeeds
        seed = eaSeeds(s);
        fprintf('  Seed %3d (%d/%d) ...', seed, s, numSeeds);
        rng(seed);

        tic;
        [result, history] = ea_lqr_codesign(A, B, ea_params, opts);
        elapsedAll(s, c) = toc;

        bestCostAll(:, s, c) = history.bestCost;
        unstableAll(:, s, c) = history.unstableCount;
        if isfield(history, 'repairCount')
            repairAll(:, s, c) = history.repairCount;
        end
        finalCostAll(s, c) = result.cost;

        % Convergence generation (5% tolerance)
        curve = history.bestCost;
        fv = curve(end);
        idx = find(curve <= fv * (1 + convergeTol), 1, 'first');
        if isempty(idx), idx = maxGen; end
        convGenAll(s, c) = idx;

        % Exponential rate fit
        res = curve - fv;
        vi = res > 0;
        if sum(vi) > 10
            gv = find(vi);
            p = polyfit(gv(:), log(res(vi)), 1);
            convRateAll(s, c) = -p(1);
        else
            convRateAll(s, c) = NaN;
        end

        fprintf(' %.1fs, cost=%.3f, conv=%d, lam=%.4f\n', ...
            elapsedAll(s,c), finalCostAll(s,c), convGenAll(s,c), convRateAll(s,c));
    end
    fprintf('\n');
end

%% ===================== Summary Table =====================
fprintf('==================== Summary (mean +/- std over %d seeds) ====================\n', numSeeds);
fprintf('%-28s  %-16s  %-16s\n', 'Metric', condNames{1}, condNames{2});
fprintf('%s\n', repmat('-', 1, 64));
fprintf('%-28s  %7.3f +/- %-5.3f  %7.3f +/- %-5.3f\n', 'Final cost', ...
    mean(finalCostAll(:,1)), std(finalCostAll(:,1)), ...
    mean(finalCostAll(:,2)), std(finalCostAll(:,2)));
fprintf('%-28s  %7.1f +/- %-5.1f  %7.1f +/- %-5.1f\n', 'Convergence gen (5%)', ...
    mean(convGenAll(:,1)), std(convGenAll(:,1)), ...
    mean(convGenAll(:,2)), std(convGenAll(:,2)));
fprintf('%-28s  %7.4f +/- %-6.4f  %7.4f +/- %-6.4f\n', 'Exp rate lambda', ...
    nanmean(convRateAll(:,1)), nanstd(convRateAll(:,1)), ...
    nanmean(convRateAll(:,2)), nanstd(convRateAll(:,2)));
fprintf('%-28s  %7.1f +/- %-5.1f  %7.1f +/- %-5.1f\n', 'Avg unstable (last 20 gen)', ...
    mean(mean(unstableAll(end-19:end,:,1),1)), std(mean(unstableAll(end-19:end,:,1),1)), ...
    mean(mean(unstableAll(end-19:end,:,2),1)), std(mean(unstableAll(end-19:end,:,2),1)));
fprintf('%-28s  %7.1f +/- %-5.1f  %7.1f +/- %-5.1f\n', 'Runtime (s)', ...
    mean(elapsedAll(:,1)), std(elapsedAll(:,1)), ...
    mean(elapsedAll(:,2)), std(elapsedAll(:,2)));
fprintf('%s\n', repmat('=', 1, 64));

%% ===================== Plotting =====================
gens = (1:maxGen)';
cb = [0.00 0.45 0.74;    % blue  (baseline)
      0.85 0.33 0.10];   % red-orange (Gershgorin)
ca = [0.75 0.85 1.0;     % light blue  (envelope)
      1.00 0.82 0.72];   % light orange (envelope)

fig = figure('Position', [60 60 1100 820], 'Color', 'w');

% --- Panel (a): Mean Best Cost ---
subplot(2, 2, 1); hold on;
for c = 1:nCond
    mu = mean(bestCostAll(:,:,c), 2);
    sd = std(bestCostAll(:,:,c), 0, 2);
    fill([gens; flipud(gens)], [mu+sd; flipud(mu-sd)], ...
        ca(c,:), 'EdgeColor', 'none', 'FaceAlpha', 0.45);
end
h = gobjects(nCond, 1);
for c = 1:nCond
    mu = mean(bestCostAll(:,:,c), 2);
    h(c) = plot(gens, mu, '-', 'Color', cb(c,:), 'LineWidth', 2);
end
xlabel('Generation'); ylabel('Best Cost');
title('(a) Best Cost (mean \pm std)');
legend(h, condNames, 'Location', 'northeast', 'FontSize', 8);
grid on; box on;

% --- Panel (b): Mean Unstable Count ---
subplot(2, 2, 2); hold on;
for c = 1:nCond
    mu = mean(unstableAll(:,:,c), 2);
    sd = std(unstableAll(:,:,c), 0, 2);
    fill([gens; flipud(gens)], [mu+sd; flipud(mu-sd)], ...
        ca(c,:), 'EdgeColor', 'none', 'FaceAlpha', 0.45);
end
for c = 1:nCond
    mu = mean(unstableAll(:,:,c), 2);
    plot(gens, mu, '-', 'Color', cb(c,:), 'LineWidth', 2);
end
xlabel('Generation'); ylabel('Unstable Individuals');
title('(b) Stability Failures (mean \pm std)');
legend(condNames, 'Location', 'northeast', 'FontSize', 8);
ylim([0 ea_params.popSize]);
grid on; box on;

% --- Panel (c): Mean Repair Count (Gershgorin only) ---
subplot(2, 2, 3); hold on;
mu_rep = mean(repairAll(:,:,2), 2);
sd_rep = std(repairAll(:,:,2), 0, 2);
fill([gens; flipud(gens)], [mu_rep+sd_rep; flipud(max(mu_rep-sd_rep,0))], ...
    ca(2,:), 'EdgeColor', 'none', 'FaceAlpha', 0.45);
plot(gens, mu_rep, '-', 'Color', cb(2,:), 'LineWidth', 2);
xlabel('Generation'); ylabel('Repairs per Generation');
title('(c) Gershgorin Repairs Applied');
grid on; box on;

% --- Panel (d): Final Cost Comparison ---
subplot(2, 2, 4); hold on;
grp = [ones(numSeeds,1); 2*ones(numSeeds,1)];
dat = [finalCostAll(:,1); finalCostAll(:,2)];
boxplot(dat, grp, 'Labels', condNames, 'Widths', 0.45, ...
    'Colors', cb, 'Symbol', '');
% overlay individual points
jitter = 0.08;
for c = 1:nCond
    x = c + jitter * randn(numSeeds, 1);
    scatter(x, finalCostAll(:,c), 40, cb(c,:), 'filled', ...
        'MarkerFaceAlpha', 0.7);
end
ylabel('Final Cost');
title('(d) Final Cost Distribution');
grid on; box on;

sgtitle(sprintf('EA-LQR Co-Design: Baseline vs Gershgorin Repair  (Grid %d\\times%d, %d seeds)', ...
    gridSize, gridSize, numSeeds), 'FontSize', 13, 'FontWeight', 'bold');

%% ===================== Save =====================
if ~exist('results', 'dir'), mkdir('results'); end

saveName = sprintf('results/gers_comparison_grid%d_%dseeds.mat', gridSize, numSeeds);
save(saveName, 'eaSeeds', 'bestCostAll', 'unstableAll', 'repairAll', ...
    'finalCostAll', 'convGenAll', 'convRateAll', 'elapsedAll', ...
    'ea_params', 'opts_base', 'gridSize', 'condNames');
fprintf('\nData saved to %s\n', saveName);

saveas(fig, sprintf('results/gers_comparison_grid%d_%dseeds.png', gridSize, numSeeds));
fprintf('Figure saved to results/\n');
fprintf('\n========== Experiment Complete ==========\n');
