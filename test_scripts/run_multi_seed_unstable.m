%% Multi-Seed A/B Experiment on an OPEN-LOOP UNSTABLE system
% Amplify spectral radius of A so rho(A)>1, then compare EA
% with vs without Gershgorin repair.  3 seeds, averaged results.
%
% Author: Pengyang Wu
% Date: 2026-02-03

clear; clc; close all;
addpath(genpath(pwd));

%% ===================== Configuration =====================
gridSize       = 5;
connectThresh  = 0.5;
Ts             = 0.2;
actDensity     = 1;
topoSeed       = 10;

eaSeeds  = [1 10];
numSeeds = length(eaSeeds);

ea_params.popSize  = 20;
ea_params.maxGen   = 150;
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

%% ===================== System Generation (unstable A) =====================
fprintf('Generating system (grid %dx%d, topoSeed=%d)...\n', gridSize, gridSize, topoSeed);
numNodes = gridSize * gridSize;
[adjMtx, nodeCoords, susceptMtx, inertiasInv, dampings] = ...
    generate_grid_topology(gridSize, connectThresh, topoSeed);
rng(topoSeed);
numActs       = round(actDensity * numNodes);
actuatedNodes = randsample(numNodes, numActs);
sys = generate_grid_plant(actuatedNodes, adjMtx, susceptMtx, inertiasInv, dampings, Ts);
A_orig = sys.A;
B = sys.B2;
[Nx, Nu] = size(B);

% --- Amplify spectral radius to make A open-loop unstable ---
rho_orig = max(abs(eig(A_orig)));
rho_target = 1.1;                       % desired spectral radius > 1
scale = rho_target / rho_orig;
A = A_orig * scale;

rho_A = max(abs(eig(A)));
fprintf('System: Nx=%d, Nu=%d\n', Nx, Nu);
fprintf('Original rho(A) = %.4f\n', rho_orig);
fprintf('Scaled   rho(A) = %.4f  (target %.2f)\n\n', rho_A, rho_target);

%% ===================== Pre-allocate =====================
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
    fprintf('========== Condition %d: %s ==========\n', c, condNames{c});

    for s = 1:numSeeds
        seed = eaSeeds(s);
        fprintf('  Seed %3d (%d/%d) ...', seed, s, numSeeds);
        rng(seed);

        tic;
        [result, history] = ea_lqr_codesign_gershgorin(A, B, ea_params, opts);
        elapsedAll(s, c) = toc;

        bestCostAll(:, s, c) = history.bestCost;
        unstableAll(:, s, c) = history.unstableCount;
        if isfield(history, 'repairCount')
            repairAll(:, s, c) = history.repairCount;
        end
        finalCostAll(s, c) = result.cost;

        % Convergence generation
        curve = history.bestCost;
        fv = curve(end);
        idx = find(curve <= fv * (1 + convergeTol), 1, 'first');
        if isempty(idx), idx = maxGen; end
        convGenAll(s, c) = idx;

        % Exponential rate
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
fprintf('============= Summary (mean +/- std, %d seeds, rho(A)=%.2f) =============\n', numSeeds, rho_A);
fprintf('%-28s  %-16s  %-16s\n', 'Metric', condNames{1}, condNames{2});
fprintf('%s\n', repmat('-', 1, 64));
printRow = @(name, v1, v2) fprintf('%-28s  %7.3f +/- %-5.3f  %7.3f +/- %-5.3f\n', ...
    name, mean(v1), std(v1), mean(v2), std(v2));
printRow('Final cost',           finalCostAll(:,1), finalCostAll(:,2));
printRow('Convergence gen (5%)', convGenAll(:,1),   convGenAll(:,2));
fprintf('%-28s  %7.4f +/- %-6.4f  %7.4f +/- %-6.4f\n', 'Exp rate lambda', ...
    nanmean(convRateAll(:,1)), nanstd(convRateAll(:,1)), ...
    nanmean(convRateAll(:,2)), nanstd(convRateAll(:,2)));
printRow('Avg unstable (last 20)', ...
    mean(unstableAll(end-19:end,:,1),1)', mean(unstableAll(end-19:end,:,2),1)');
printRow('Runtime (s)', elapsedAll(:,1), elapsedAll(:,2));
fprintf('%s\n', repmat('=', 1, 64));

%% ===================== Dense LQR Reference (for normalisation) =====================
Q_bm_ref = eye(Nx);  R_bm_ref = eye(Nu);
K_dense_ref = -dlqr(A, B, Q_bm_ref, R_bm_ref);
K_dense_ref(abs(K_dense_ref) < 1e-3) = 0;
costBM_ref  = get_lqr_cost(A, B, Q_bm_ref, R_bm_ref, K_dense_ref);
alpha_ref   = ea_params.alpha;
w_c_ref = 0.05*(1-alpha_ref);  w_r_ref = 0.4*(1-alpha_ref);  w_s_ref = 0.2*(1-alpha_ref);
J_dense_ref = 1.0 + w_c_ref*nnz(K_dense_ref) ...
                  + w_r_ref*nnz(any(K_dense_ref,2)) ...
                  + w_s_ref*nnz(any(K_dense_ref,1));
fprintf('Dense LQR ref: costBM=%.4e  J_dense_total=%.4f\n', costBM_ref, J_dense_ref);

%% ===================== Figure =====================
gens = (1:maxGen)';
cb = [0.00 0.45 0.74;    % blue  (baseline)
      0.85 0.33 0.10];   % orange (Gershgorin)
ca = [0.75 0.85 1.0;
      1.00 0.82 0.72];

% Font sizes for double-column LaTeX
axFS  = 13;   % axis tick labels
labFS = 14;   % axis labels
titFS = 15;   % subplot titles
legFS = 11;   % legend
lw    = 2.5;  % line width

fig = figure('Position', [60 60 1400 420], 'Color', 'w');
set(fig, 'DefaultAxesFontSize', axFS, 'DefaultAxesFontName', 'Times New Roman');

% --- (a) Best Cost (normalised by Dense LQR total cost) ---
subplot(1, 3, 1); hold on;
for c = 1:nCond
    mu = mean(bestCostAll(:,:,c), 2) / J_dense_ref;
    sd = std(bestCostAll(:,:,c), 0, 2)  / J_dense_ref;
    fill([gens; flipud(gens)], [mu+sd; flipud(mu-sd)], ...
        ca(c,:), 'EdgeColor', 'none', 'FaceAlpha', 0.45);
end
h = gobjects(nCond,1);
for c = 1:nCond
    h(c) = plot(gens, mean(bestCostAll(:,:,c),2) / J_dense_ref, '-', ...
        'Color', cb(c,:), 'LineWidth', lw);
end
hRef = yline(1.0, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.8);
xlabel('Generation', 'FontSize', labFS);
ylabel('Best Cost $/ \; J_{\mathrm{dense}}$', 'FontSize', labFS, 'Interpreter', 'latex');
title('(a) Best Cost Convergence', 'FontSize', titFS);
legend([h; hRef], [condNames, {'Dense LQR ($=1$)'}], ...
    'Location', 'northeast', 'FontSize', legFS, 'Interpreter', 'latex');
grid on; box on;

% --- (b) Unstable Count ---
subplot(1, 3, 2); hold on;
for c = 1:nCond
    mu = mean(unstableAll(:,:,c), 2);
    sd = std(unstableAll(:,:,c), 0, 2);
    fill([gens; flipud(gens)], [mu+sd; flipud(mu-sd)], ...
        ca(c,:), 'EdgeColor', 'none', 'FaceAlpha', 0.45);
end
for c = 1:nCond
    plot(gens, mean(unstableAll(:,:,c),2), '-', 'Color', cb(c,:), 'LineWidth', lw);
end
xlabel('Generation', 'FontSize', labFS);
ylabel('Unstable Count', 'FontSize', labFS);
title('(b) Stability Failures', 'FontSize', titFS);
legend(condNames, 'Location', 'northeast', 'FontSize', legFS);
ylim([0 ea_params.popSize]);
grid on; box on;

% --- (c) Gershgorin Repairs ---
subplot(1, 3, 3); hold on;
mu_r = mean(repairAll(:,:,2), 2);
sd_r = std(repairAll(:,:,2), 0, 2);
fill([gens; flipud(gens)], [mu_r+sd_r; flipud(max(mu_r-sd_r,0))], ...
    ca(2,:), 'EdgeColor', 'none', 'FaceAlpha', 0.45);
plot(gens, mu_r, '-', 'Color', cb(2,:), 'LineWidth', lw);
xlabel('Generation', 'FontSize', labFS);
ylabel('Repairs / Gen', 'FontSize', labFS);
title('(c) Gershgorin Repairs', 'FontSize', titFS);
grid on; box on;

sgtitle(sprintf(['EA-LQR: Baseline vs Gershgorin  ' ...
    '(Grid %d\\times%d, \\rho(A)=%.2f, %d seeds)'], ...
    gridSize, gridSize, rho_A, numSeeds), ...
    'FontSize', 16, 'FontWeight', 'bold');

%% ===================== Save =====================
if ~exist('results', 'dir'), mkdir('results'); end
tag = sprintf('grid%d_unstable_%dseeds', gridSize, numSeeds);
save(fullfile('results', ['gers_comparison_' tag '.mat']), ...
    'eaSeeds', 'bestCostAll', 'unstableAll', 'repairAll', ...
    'finalCostAll', 'convGenAll', 'convRateAll', 'elapsedAll', ...
    'ea_params', 'opts_base', 'gridSize', 'rho_A', 'condNames');
saveas(fig, fullfile('results', ['gers_comparison_' tag '.png']));
fprintf('\nResults saved to results/\n');
fprintf('========== Experiment Complete ==========\n');

%% ---- Save plot-only data ----
plotDir = fullfile(pwd, 'plot_only');
if ~exist(plotDir, 'dir'), mkdir(plotDir); end
save(fullfile(plotDir, 'fig_data_unstable.mat'), ...
    'bestCostAll', 'unstableAll', 'repairAll', ...
    'J_dense_ref', 'ea_params', 'maxGen', 'nCond', 'condNames', ...
    'gridSize', 'rho_A', 'numSeeds', 'eaSeeds');
fprintf('Unstable plot data saved → %s/fig_data_unstable.mat\n', plotDir);
