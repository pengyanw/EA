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

eaSeeds  = [1 10 15];
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

% --- Which variant of Algorithm 1 to run (set explicitly, do not inherit) ---
% Section III of the paper: every individual starts at ell_0 = ||K_d||_0, so the
% cost curve begins at the dense controller (normalized ~1) and prunes downward.
% Set initLinks='random' to recover the earlier behaviour, in which ell_0 was
% drawn uniformly and the curve started from an arbitrary sparsity instead.
% Keep this consistent with analysis_perf_bounds.m (Figure 1).
opts_base.initLinks        = 'paper';    % 'random' => ell_0 ~ U{1, ||K_d||_0}
opts_base.gateLinkMutation = false;      % true     => delta gated on p_m

% Canonicalise the gene's link count after each decode (see analysis_perf_bounds.m
% and ea_lqr_codesign_gershgorin.m). NOTE the interaction with repair: when the
% Gershgorin branch succeeds it overwrites ell with nnz(K_sparse) anyway, so on
% repaired individuals that write-back wins over the canonicalisation. The two
% rules disagree -- nnz counts surviving entries, the canonical ell is the
% deepest surviving rank -- and only the latter leaves K unchanged.
opts_base.linkDecode       = 'canonical';

% Draw the greedy baseline (Algorithm 3, EA functions/greedy_prune.m) in panel
% (a), one curve per repair condition. Set false to omit both curves and skip
% the runs that produce them. Greedy is deterministic, so one run per condition
% suffices and there is no shaded band. This is the figure the "repair separates
% the two searches" claim rests on: repair improves the EA and degrades greedy,
% which is not visible from the EA curves alone.
showGreedy = false;

% Repair step size: Inf is the pure Polyak step of Eq. (37), which is what
% Proposition 5 proves. It is NOT equivalent to the previous clamped step at
% full scale -- measured over 15,415 repair calls it lowers the repair failure
% rate from 42.0% to 38.1% (final cost 21.49 -> 21.38, within one sigma).
opts_base.gersEtaCap = Inf;

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
rho_target = 1.3;                       % desired spectral radius > 1
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
phaseTot     = zeros(5, nCond);   % rows: repair exit phase -1, 0, 1, 2, 3

%% ===================== Run Experiments =====================
gdCurveCond = NaN(maxGen, nCond);   % greedy incumbent on the EA generation axis
gdFinalCond = NaN(1, nCond);

for c = 1:nCond
    opts = opts_base;
    opts.useGersRepair = (c == 2);
    fprintf('========== Condition %d: %s ==========\n', c, condNames{c});

    if showGreedy
        % Same gene, same decode, same cost oracle, same repair setting.
        [gdB, gdH] = greedy_prune(A, B, ea_params, ...
            struct('Q', eye(Nx), 'R', eye(Nu), 'denseTol', 1e-3, ...
                   'useGersRepair', opts.useGersRepair));
        gdCurveCond(:, c) = gdH.JperGen;
        gdFinalCond(c)    = gdB.J;
        fprintf('  greedy: J=%.4f  N_a=%d N_c=%d  %d evals (%.2fx EA budget), %d infeasible\n', ...
            gdB.J, gdB.Na, gdB.Nc, gdB.evals, ...
            gdB.evals/(ea_params.popSize*maxGen), gdB.nInf);
    end

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
        % Exit-phase mix of Algorithm 2, reported in Section VII:
        % rows are phase -1 (failed) / 0 / 1 / 2 / 3.
        if isfield(history, 'repairPhase')
            phaseTot(:, c) = phaseTot(:, c) + sum(history.repairPhase, 2);
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

%% ===================== Repair exit-phase mix (for Section VII) =====================
tot = sum(phaseTot(:,2));
if tot > 0
    lbl = {'failed', 'phase 0 (already stable)', ...
           'phase 1 (subgradient, in K_S)', ...
           'phase 2 (line search, in K_S)', ...
           'phase 3 (support re-selection)'};
    fprintf('\nAlgorithm 2 exit phases over %d repair calls:\n', tot);
    for k = 1:5
        fprintf('  %-32s %7d  (%5.2f%%)\n', lbl{k}, phaseTot(k,2), 100*phaseTot(k,2)/tot);
    end
    ok = tot - phaseTot(1,2);
    if ok > 0
        fprintf('  of the %d SUCCESSFUL repairs: phase 1 %.1f%%, phase 2 %.1f%%, phase 3 %.1f%%\n', ...
            ok, 100*phaseTot(3,2)/ok, 100*phaseTot(4,2)/ok, 100*phaseTot(5,2)/ok);
    end
end

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
axFS  = 17;   % axis tick labels
labFS = 19;   % axis labels
titFS = 19;   % subplot titles
legFS = 15;   % legend
lw    = 2.5;  % line width

% Two panels, tiled with compact spacing so the exported PDF has little dead
% space. Panel (c) (repairs per generation) was removed; repairAll is still
% saved to the .mat and plot_only/plot_fig2_pdf.m can rebuild any variant.
fig = figure('Position', [60 60 1040 400], 'Color', 'w');
% Times for tick labels, and DefaultTextFontName so that titles, axis labels and
% legend entries pick it up too -- MATLAB's 'latex' interpreter always renders
% Computer Modern, so every string below stays on the 'tex' interpreter.
set(fig, 'DefaultAxesFontSize', axFS, 'DefaultAxesFontName', 'Times New Roman', ...
         'DefaultTextFontName', 'Times New Roman', ...
         'DefaultLegendFontName', 'Times New Roman');
tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'tight', 'Padding', 'tight');

% --- (a) Best Cost (normalised by Dense LQR total cost) ---
nexttile(tl); hold on;
for c = 1:nCond
    mu = mean(bestCostAll(:,:,c), 2) / J_dense_ref;
    sd = std(bestCostAll(:,:,c), 0, 2)  / J_dense_ref;
    fill([gens; flipud(gens)], [mu+sd; flipud(mu-sd)], ...
        ca(c,:), 'EdgeColor', 'none', 'FaceAlpha', 0.45);
end
ha = gobjects(nCond,1);  hg = gobjects(nCond,1);
for c = 1:nCond
    ha(c) = plot(gens, mean(bestCostAll(:,:,c),2) / J_dense_ref, '-', ...
        'Color', cb(c,:), 'LineWidth', lw);
end
if showGreedy
    for c = 1:nCond
        hg(c) = plot(gens, gdCurveCond(:,c) / J_dense_ref, ':', ...
            'Color', cb(c,:), 'LineWidth', lw*0.9);
    end
    legend([ha; hg], [strcat("EA, ", condNames), strcat("Greedy, ", condNames)], ...
        'FontSize', legFS, 'Location', 'northeast', 'Box', 'off', ...
        'AutoUpdate', 'off');
end
% HandleVisibility off, or the dense reference is auto-appended to the legend
% above as an unnamed entry.
yline(1.0, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.8, 'HandleVisibility', 'off');
xlabel('Generation', 'FontSize', labFS);
ylabel('Best cost / {\itJ}_{dense}', 'FontSize', labFS);
title('(a) Best Cost Convergence', 'FontSize', titFS);
% No legend here: panel (b) carries the same two conditions in the same colors.
% The dashed grey line at 1 is the dense-LQR reference -- state that in the
% caption, since it no longer has a legend entry.
grid on; box on;

% --- (b) Unstable Count ---
nexttile(tl); hold on;
for c = 1:nCond
    mu = mean(unstableAll(:,:,c), 2);
    sd = std(unstableAll(:,:,c), 0, 2);
    fill([gens; flipud(gens)], [mu+sd; flipud(max(mu-sd,0))], ...
        ca(c,:), 'EdgeColor', 'none', 'FaceAlpha', 0.45);
end
hb = gobjects(nCond,1);
for c = 1:nCond
    hb(c) = plot(gens, mean(unstableAll(:,:,c),2), '-', 'Color', cb(c,:), 'LineWidth', lw);
end
xlabel('Generation', 'FontSize', labFS);
ylabel('Unstable Count', 'FontSize', labFS);
title('(b) Stability Failures', 'FontSize', titFS);
% Legend on the line handles; the fill() patches are drawn first and would
% otherwise be what gets labelled.
legend(hb, condNames, 'Location', 'northeast', 'FontSize', legFS);
ylim([0 ea_params.popSize]);
grid on; box on;

% The layout title duplicates the LaTeX caption and costs a band of vertical
% space, so it is off for the paper export. Set true when the PNG is viewed on
% its own and needs to carry its own context.
standaloneTitle = false;
if standaloneTitle
    title(tl, sprintf(['EA-LQR: Baseline vs Gershgorin  ' ...
        '(Grid %d\\times%d, \\rho(A)=%.2f, %d seeds)'], ...
        gridSize, gridSize, rho_A, numSeeds), ...
        'FontSize', titFS, 'FontWeight', 'bold');
end

%% ===================== Save =====================
if ~exist('results', 'dir'), mkdir('results'); end
tag = sprintf('grid%d_unstable_%dseeds', gridSize, numSeeds);
save(fullfile('results', ['gers_comparison_' tag '.mat']), ...
    'eaSeeds', 'bestCostAll', 'unstableAll', 'repairAll', ...
    'finalCostAll', 'convGenAll', 'convRateAll', 'elapsedAll', ...
    'ea_params', 'opts_base', 'gridSize', 'rho_A', 'condNames');
% Vector PDF for the paper (Figure 2). exportgraphics with ContentType 'vector'
% keeps the text selectable and the lines resolution-independent;
% saveas(...,'pdf') would go through the printer path and pad the page.
% Exporting the tiledlayout rather than the figure crops to the drawn content.
pdfPath = fullfile('results', ['gers_comparison_' tag '.pdf']);
exportgraphics(tl, pdfPath, 'ContentType', 'vector', 'BackgroundColor', 'white');
exportgraphics(tl, fullfile('results', ['gers_comparison_' tag '.png']), ...
    'Resolution', 200);
fprintf('\nFigure 2 (vector PDF) saved to %s\n', pdfPath);
fprintf('Results saved to results/\n');
fprintf('========== Experiment Complete ==========\n');

%% ---- Save plot-only data ----
plotDir = fullfile(pwd, 'plot_only');
if ~exist(plotDir, 'dir'), mkdir(plotDir); end
% gdCurveCond/gdFinalCond are all-NaN when showGreedy is false; plot_fig2_pdf.m
% must guard with any(isfinite(...)) rather than exist().
save(fullfile(plotDir, 'fig_data_unstable.mat'), ...
    'bestCostAll', 'unstableAll', 'repairAll', ...
    'J_dense_ref', 'ea_params', 'maxGen', 'nCond', 'condNames', ...
    'gridSize', 'rho_A', 'numSeeds', 'eaSeeds', ...
    'gdCurveCond', 'gdFinalCond', 'showGreedy');
fprintf('Unstable plot data saved → %s/fig_data_unstable.mat\n', plotDir);
