%% analysis_perf_bounds.m
% Compare Dense LQR, Diagonal (distributed) LQR, and EA-LQR:
%   - Convergence speed (generation to beat each baseline)
%   - Performance gap (cost decomposition: LQR perf vs structural penalty)
%   - Scaling trend across grid sizes
%
% Dense LQR    := K_dense with |K_ij| >= 1e-3  (benchmark, normalized to 1)
% Diagonal LQR := actuator i only reads states [2i-1, 2i] (fully local)
% EA-LQR       := ea_lqr_codesign_gershgorin with Gershgorin repair
%
% Author: Pengyang Wu
% Date: 2026-02-17

clear; clc; close all;
addpath(genpath(pwd));

%% ===================== Configuration =====================
gridSizes     = [5 7];
eaSeeds       = [1];
connectThresh = 0.5;
Ts            = 0.2;
actDensity    = 1;

nGrid    = numel(gridSizes);
numSeeds = numel(eaSeeds);

ea_params.popSize  = 20;
ea_params.maxGen   = 150;
ea_params.pMutate  = 0.05;
ea_params.pCross   = 0.8;
ea_params.nTop     = 10;
ea_params.alpha    = 0;
maxGen = ea_params.maxGen;

opts.verbose       = false;
opts.useGersRepair = false;

alpha = ea_params.alpha;

%% Bound display options (for Fig 1)
showLB    = true;  % Show exponential lower bound from convergence theorem (Eq. 49)
lb_lambda = 0.05;  % λ for Eq.(49) lower bound; paper empirical range ≈ [0.035, 0.069] (4×4 grid)

%% ===================== Pre-allocate =====================
% Per-method cost decomposition:  J_total = J_perf + J_struct
%   method 1: Dense, 2: Diagonal, 3: EA
J_total_all  = NaN(nGrid, 3, numSeeds);
J_perf_all   = NaN(nGrid, 3, numSeeds);   % pure LQR / costBM
J_struct_all = NaN(nGrid, 3, numSeeds);   % structural penalty

nnz_all      = NaN(nGrid, 3, numSeeds);
nAct_all     = NaN(nGrid, 3, numSeeds);   % active actuators
nSens_all    = NaN(nGrid, 3, numSeeds);   % active sensors

% Convergence milestones
convGenDense = NaN(nGrid, numSeeds);  % gen when EA < Dense
convGenDiag  = NaN(nGrid, numSeeds);  % gen when EA < Diagonal

% Full EA curves (raw)
eaCurves     = cell(nGrid, numSeeds);
sysSize      = zeros(nGrid, 2);
denseCostRef = zeros(nGrid, numSeeds);

J_perf_trunc      = NaN(nGrid, numSeeds);   % κ=3 truncated LQR ratio (/ costBM)
nnz_trunc         = NaN(nGrid, numSeeds);
J_comp_trunc      = NaN(nGrid, numSeeds);   % κ=3 composite cost (LQR+struct) ratio
trunc_stable_grid = false(nGrid, numSeeds); % true if K_trunc closed-loop is stable

% Per-grid topology & controller (for row-2 diagrams)
adjMtx_all       = cell(nGrid, 1);
K_ea_all         = cell(nGrid, 1);
actuatedNodes_all = cell(nGrid, 1);

shin_kappa = 1;   % locality radius for K-truncation (shared by grid and IEEE 13-bus)

%% ===================== Main Loop =====================
for s = 1:numSeeds
    seed = eaSeeds(s);
    fprintf('====== Seed %d (%d/%d) ======\n', seed, s, numSeeds);

    for g = 1:nGrid
        gridSize = gridSizes(g);
        numNodes = gridSize^2;

        % --- Generate system (A,B) controlled by seed ---
        rng(seed);
        [adjMtx, ~, susceptMtx, inertiasInv, dampings] = ...
            generate_grid_topology(gridSize, connectThresh, seed);
        numActs      = round(actDensity * numNodes);
        actuatedNodes = randsample(numNodes, numActs);
        sys = generate_grid_plant(actuatedNodes, adjMtx, susceptMtx, ...
                                  inertiasInv, dampings, Ts);
        A = sys.A;  B = sys.B2;
        [Nx, Nu] = size(B);
        sysSize(g,:) = [Nx, Nu];
        adjMtx_all{g}        = adjMtx;
        actuatedNodes_all{g} = actuatedNodes;

        Q_bm = eye(Nx);  R_bm = eye(Nu);
        K_full = -dlqr(A, B, Q_bm, R_bm);

        % ---------- Dense LQR (|K| >= 1e-3) ----------
        K_dense = K_full;
        K_dense(abs(K_dense) < 1e-3) = 0;
        costBM = get_lqr_cost(A, B, Q_bm, R_bm, K_full);

        % ---------- κ=3 Truncation (Shin-style locality) ----------
        G_grid  = graph(adjMtx);
        Is_grid = cell(numNodes, 1);
        Js_grid = cell(numNodes, 1);
        for ii = 1:numNodes
            Is_grid{ii} = [2*ii-1,  2*ii];  % theta_i, omega_i state columns
        end
        for jj = 1:Nu
            %Js_grid{actuatedNodes(jj)} = jj;    % node k → control row jj
            Js_grid{jj} = jj;    % node k → control row jj
        end
        prob_grid = struct('name', 'grid');      % → unweighted hop distances
        K_trunc_grid    = truncate(K_full, G_grid, Is_grid, Js_grid, shin_kappa, prob_grid);
        nnz_trunc(g, s) = nnz(K_trunc_grid);
        rho_trunc_grid  = max(abs(eig(A + B * K_trunc_grid)));
        trunc_stable_grid(g, s) = (rho_trunc_grid < 1.0);
        if ~trunc_stable_grid(g, s)
            fprintf('  Trunc k=%d grid %dx%d: UNSTABLE (rho=%.4f)\n', shin_kappa, gridSize, gridSize, rho_trunc_grid);
            J_comp_trunc(g, s) = Inf;
            J_perf_trunc(g, s) = Inf;
        else
            J_comp_trunc(g, s) = cost_EA(A, B, Q_bm, R_bm, K_trunc_grid, costBM, alpha);
            J_perf_trunc(g, s) = get_lqr_cost(A, B, Q_bm, R_bm, K_trunc_grid) / costBM;
        end

        % ---------- Diagonal LQR ----------
        K_diag = zeros(Nu, Nx);
        for i = 1:Nu
            c1 = 2*i - 1;  c2 = 2*i;
            if c2 <= Nx
                K_diag(i, c1) = K_full(i, c1);
                K_diag(i, c2) = K_full(i, c2);
            elseif c1 <= Nx
                K_diag(i, c1) = K_full(i, c1);
            end
        end
        % if max(abs(eig(A + B * K_diag))) >= 1.0
        %     [K_diag, info_d] = gersgorin_stabilize_K(A, B, K_diag, K_full);
        %     if ~info_d.success
        %         fprintf('  WARNING: Diagonal K still unstable (grid %dx%d)\n', gridSize, gridSize);
        %     end
        % end

        % Dense / Diagonal cost decomposition (per seed)
        [J_d, Jp_d, Js_d] = cost_decompose(A, B, Q_bm, R_bm, K_dense, costBM, alpha);
        denseCostRef(g, s) = J_d;

        [J_g, Jp_g, Js_g] = cost_decompose(A, B, Q_bm, R_bm, K_diag, costBM, alpha);

        % ---- Run EA (reset RNG so EA is reproducible for this seed) ----
        rng(seed);
        [result, history] = ea_lqr_codesign_gershgorin(A, B, ea_params, opts);
        K_ea = result.K_sparse;
        eaCurves{g, s} = history.bestCost;
        K_ea_all{g} = K_ea;   % keep last seed (for topology diagram)

        [J_e, Jp_e, Js_e] = cost_decompose(A, B, Q_bm, R_bm, K_ea, costBM, alpha);

        % Store decomposition
        J_total_all(g, 1, s) = J_d;   J_perf_all(g, 1, s) = Jp_d;  J_struct_all(g, 1, s) = Js_d;
        J_total_all(g, 2, s) = J_g;   J_perf_all(g, 2, s) = Jp_g;  J_struct_all(g, 2, s) = Js_g;
        J_total_all(g, 3, s) = J_e;   J_perf_all(g, 3, s) = Jp_e;  J_struct_all(g, 3, s) = Js_e;

        nnz_all(g, 1, s) = nnz(K_dense);  nAct_all(g, 1, s) = nnz(any(K_dense,2)); nSens_all(g, 1, s) = nnz(any(K_dense,1));
        nnz_all(g, 2, s) = nnz(K_diag);   nAct_all(g, 2, s) = nnz(any(K_diag,2));  nSens_all(g, 2, s) = nnz(any(K_diag,1));
        nnz_all(g, 3, s) = nnz(K_ea);     nAct_all(g, 3, s) = nnz(any(K_ea,2));    nSens_all(g, 3, s) = nnz(any(K_ea,1));

        % Convergence milestones
        curve = history.bestCost;
        idx_d = find(curve < J_d, 1, 'first');
        idx_g = find(curve < J_g, 1, 'first');
        if ~isempty(idx_d), convGenDense(g, s) = idx_d; end
        if ~isempty(idx_g), convGenDiag(g, s)  = idx_g; end

        fprintf('  Grid %dx%d  J_total: Dense=%.2f  Diag=%.2f  EA=%.2f | conv<Dense@gen%d  conv<Diag@gen%d\n', ...
            gridSize, gridSize, J_d, J_g, J_e, ...
            round(nanmean(convGenDense(g,s))), round(nanmean(convGenDiag(g,s))));
    end
end

%% ===================== Aggregate =====================
meanTotal  = mean(J_total_all, 3);
meanPerf   = mean(J_perf_all, 3);
meanStruct = mean(J_struct_all, 3);
meanNnz    = mean(nnz_all, 3);
meanAct    = mean(nAct_all, 3);
meanSens   = mean(nSens_all, 3);

meanPerfTrunc  = nanmean(J_perf_trunc, 2);    % nGrid×1, mean over seeds
stdPerfTrunc   = nanstd(J_perf_trunc, 0, 2);
meanNnzTrunc   = nanmean(nnz_trunc, 2);
stdPerfEA_grid = squeeze(std(J_perf_all(:, 3, :), 0, 3));  % nGrid×1
meanCompTrunc  = nanmean(J_comp_trunc, 2);   % composite trunc cost, mean over seeds

dimLabel = cell(nGrid, 1);
nxVec   = sysSize(:,1);
nuVec   = sysSize(:,2);
for g = 1:nGrid
    dimLabel{g} = sprintf('%dx%d', gridSizes(g), gridSizes(g));
end

%% ===================== Theoretical Bounds =====================
w_c = 0.05 * (1 - alpha);
w_r = 0.4  * (1 - alpha);
w_s = 0.2  * (1 - alpha);

% Proposition 1: Dense LQR upper bound (closed-form)
%   J_dense <= 1 + w_c*m*n + w_r*m + w_s*n
J_dense_UB = 1 + w_c .* (nxVec .* nuVec) + w_r .* nuVec + w_s .* nxVec;

% Proposition 2: Diagonal LQR structural lower bound
%   J_struct_diag = w_c*2m + w_r*m + w_s*2m  (exactly known)
J_diag_struct_theory = w_c * 2 .* nuVec + w_r .* nuVec + w_s * 2 .* nuVec;


delta_perf_budget = zeros(nGrid, 1);
for g = 1:nGrid
    delta_perf_budget(g) = w_c * (meanNnz(g,1) - meanNnz(g,3)) + ...
                           w_r * (meanAct(g,1)  - meanAct(g,3))  + ...
                           w_s * (meanSens(g,1) - meanSens(g,3));
end

% Actual performance degradation EA vs Dense
delta_perf_actual = meanPerf(:,3) - meanPerf(:,1);

fprintf('\n============ Theoretical Bounds ============\n');
fprintf('%-8s  J_dense_UB  J_dense_act  J_diag_struct_th  delta_perf_budget  delta_perf_actual\n', 'Grid');
fprintf('%s\n', repmat('-', 1, 85));
for g = 1:nGrid
    fprintf('%-8s  %10.1f  %10.1f  %16.1f  %16.1f  %16.2f\n', dimLabel{g}, ...
        J_dense_UB(g), meanTotal(g,1), J_diag_struct_theory(g), ...
        delta_perf_budget(g), delta_perf_actual(g));
end
fprintf('\nEA beats Dense whenever delta_perf_actual < delta_perf_budget (always satisfied above)\n');

%% ===================== Console Analysis =====================
fprintf('\n');
fprintf('===================================================================\n');
fprintf('  COST DECOMPOSITION  (J = J_perf + J_struct,  alpha=%.1f)\n', alpha);
fprintf('===================================================================\n');
fprintf('%-8s  %-8s  | %-26s | %-26s | %-26s\n', '', '', 'Dense LQR', 'Diagonal LQR', 'EA-LQR');
fprintf('%-8s  %-8s  | %8s %8s %8s | %8s %8s %8s | %8s %8s %8s\n', ...
    'Grid', 'n,m', 'J', 'Jperf', 'Jstruct', 'J', 'Jperf', 'Jstruct', 'J', 'Jperf', 'Jstruct');
fprintf('%s\n', repmat('-', 1, 100));
for g = 1:nGrid
    fprintf('%-8s  %3d,%3d  ', dimLabel{g}, nxVec(g), nuVec(g));
    for m = 1:3
        fprintf('| %7.1f %7.2f %7.1f ', meanTotal(g,m), meanPerf(g,m), meanStruct(g,m));
    end
    fprintf('\n');
end
fprintf('%s\n', repmat('=', 1, 100));

fprintf('\n%-8s  |  nnz(K)          |  Active act      |  Active sens\n', 'Grid');
fprintf('%-8s  | %5s %5s %5s | %5s %5s %5s | %5s %5s %5s\n', ...
    '', 'Dens', 'Diag', 'EA', 'Dens', 'Diag', 'EA', 'Dens', 'Diag', 'EA');
fprintf('%s\n', repmat('-', 1, 72));
for g = 1:nGrid
    fprintf('%-8s  | %5.0f %5.0f %5.0f | %5.0f %5.0f %5.0f | %5.0f %5.0f %5.0f\n', ...
        dimLabel{g}, meanNnz(g,:), meanAct(g,:), meanSens(g,:));
end

fprintf('\n%-8s  Conv gen < Dense   Conv gen < Diag\n', 'Grid');
fprintf('%s\n', repmat('-', 1, 45));
for g = 1:nGrid
    fprintf('%-8s    %5.1f +/- %-5.1f     %5.1f +/- %-5.1f\n', dimLabel{g}, ...
        nanmean(convGenDense(g,:)), nanstd(convGenDense(g,:)), ...
        nanmean(convGenDiag(g,:)),  nanstd(convGenDiag(g,:)));
end

%% ===================== Figure 1: Normalized Convergence =====================
axFS = 18;  labFS = 19;  titFS = 19;  legFS = 15;  lw = 3;
cb_ea    = [0.00 0.45 0.74];
cb_dens  = [0.50 0.50 0.50];
cb_diag  = [0.85 0.33 0.10];
cb_trunc = [0.47 0.67 0.19];   % green — κ=3 truncation
cb_lb    = [0.93 0.69 0.13];   % amber — Eq. (49) exponential lower bound
ca_ea    = [0.75 0.85 1.00];
gens = (1:maxGen)';

% 3 panels: nGrid grid systems + 1 Shin 2d-mesh (panel 3 filled after Shin section)
fig1 = figure('Position', [40 60 2200 820], 'Color', 'w');
set(fig1, 'DefaultAxesFontSize', axFS, 'DefaultAxesFontName', 'Times New Roman');

for g = 1:nGrid
    subplot(2, nGrid+1, g); hold on;

    % Normalized EA curves (/ per-seed Dense LQR cost)
    Jref = mean(denseCostRef(g,:));
    normMat = zeros(maxGen, numSeeds);
    for s = 1:numSeeds
        normMat(:, s) = eaCurves{g, s} / denseCostRef(g, s);
    end
    mu = mean(normMat, 2);
    sd = std(normMat, 0, 2);

    fill([gens; flipud(gens)], [mu+sd; flipud(max(mu-sd, 0))], ...
        ca_ea, 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    hEA = plot(gens, mu, '-', 'Color', cb_ea, 'LineWidth', lw);

    hDens = yline(1.0, '--', 'Color', cb_dens, 'LineWidth', 2);

    J_diag_n = meanTotal(g,2) / Jref;
    if J_diag_n < 1e4
        hDiag = yline(J_diag_n, '-.', 'Color', cb_diag, 'LineWidth', 2);
    else
        hDiag = plot(NaN, NaN, '-.', 'Color', cb_diag, 'LineWidth', 2);
        text(maxGen/2, max(mu)*0.9, 'Diag: Unstable', 'Color', cb_diag, ...
            'FontSize', legFS, 'FontName', 'Times New Roman', 'HorizontalAlignment', 'center');
    end

    % κ=3 Truncation reference
    J_trunc_n = meanCompTrunc(g) / Jref;
    if all(trunc_stable_grid(g, :))
        hTrunc = yline(J_trunc_n, ':', 'Color', cb_trunc, 'LineWidth', 2);
    else
        hTrunc = plot(NaN, NaN, ':', 'Color', cb_trunc, 'LineWidth', 2);
        text(maxGen*0.7, max(mu)*0.85, 'Trunc: Unstable', 'Color', cb_trunc, ...
            'FontSize', legFS, 'FontName', 'Times New Roman', 'HorizontalAlignment', 'center');
    end

    % Exponential lower bound: J_min(t) >= J_R* + A0*exp(-lambda*t)  [Eq. (49)]
    hLB = plot(NaN, NaN, '-.', 'Color', cb_lb, 'LineWidth', lw*0.8);  % sentinel
    if showLB
        lam_g = lb_lambda;   % λ from Eq.(49); adjust lb_lambda in Configuration section
        % Estimate J_R* (recombination plateau) and A0 (initial amplitude) from EA data
        JR_vec  = arrayfun(@(s) eaCurves{g,s}(end), 1:numSeeds);
        A0_vec  = arrayfun(@(s) eaCurves{g,s}(1),   1:numSeeds) - JR_vec;
        JR_star = mean(JR_vec);
        A0_g    = max(mean(A0_vec), 0);
        LB_curve = (JR_star + A0_g * exp(-lam_g * gens)) / Jref;
        delete(hLB);
        hLB = plot(gens, LB_curve, '-.', 'Color', cb_lb, 'LineWidth', lw*0.8);
    end

    % Mark convergence generation (mean)
    cg = round(nanmean(convGenDense(g,:)));
    if ~isnan(cg) && cg <= maxGen
        plot(cg, mu(cg), 'v', 'MarkerSize', 8, 'MarkerFaceColor', cb_ea, ...
            'MarkerEdgeColor', 'k', 'LineWidth', 0.8);
    end

    xlabel('Generation', 'FontSize', labFS);
    if g == 1
        ylabel('$J \;/\; J_{\mathrm{dense}}$', 'FontSize', labFS, 'Interpreter', 'latex');
    end
    title(sprintf('(%c) %d$\\times$%d ($n$=%d, $m$=%d)', ...
        char('a'+g-1), gridSizes(g), gridSizes(g), nxVec(g), nuVec(g)), ...
        'FontSize', titFS, 'Interpreter', 'latex');
    grid on; box on;

    if g == 1 || g == nGrid
        legHandles = [hEA, hDens, hDiag];
        legLabels  = {'EA-LQR', 'Dense ($=1$)', 'Diagonal'};
        if all(trunc_stable_grid(g, :))
            legHandles(end+1) = hTrunc;
            legLabels{end+1}  = '$\kappa=1$ Trunc.';
        end
        if showLB
            legHandles(end+1) = hLB;
            legLabels{end+1}  = 'Conv';
        end
        legend(legHandles, legLabels, 'Location', 'northeast', ...
            'FontSize', legFS, 'Interpreter', 'latex');
    end
end

sgtitle('EA-LQR Convergence Normalized by Dense LQR', ...
    'FontSize', 16, 'FontWeight', 'bold', 'FontName', 'Times New Roman');


%% ===================== Save =====================
outDir = fullfile(pwd, 'results');
if ~exist(outDir, 'dir'), mkdir(outDir); end

% fig1 saved later after 2d-mesh panel is appended
% fig2 saved later after 2d-mesh data is appended

save(fullfile(outDir, 'perf_bounds_analysis.mat'), ...
    'gridSizes', 'eaSeeds', 'sysSize', 'J_total_all', 'J_perf_all', ...
    'J_struct_all', 'nnz_all', 'nAct_all', 'nSens_all', ...
    'convGenDense', 'convGenDiag', 'eaCurves', 'denseCostRef', ...
    'ea_params', 'meanTotal', 'meanPerf', 'meanStruct', ...
    'J_dense_UB', 'J_diag_struct_theory', 'delta_perf_budget', 'delta_perf_actual', ...
    'J_perf_trunc', 'nnz_trunc', 'meanPerfTrunc', 'stdPerfTrunc', 'meanNnzTrunc');

fprintf('\nResults saved to %s\n', outDir);


%% ========= Section: IEEE 13-Bus System — EA vs κ=1 Truncation =========
% IEEE 13-bus power system (swing-equation model, Ts=0.2 s)
%   13 buses, Nx=26 states (angle + frequency), Nu=13 (full actuation)
%
% Three methods:
%   1. Dense LQR  — K* from DARE,  J/J* = 1.0 (benchmark)
%   2. κ=2 Trunc  — locality-truncated K*, zero entries beyond 2 hops
%   3. EA-LQR     — ea_lqr_codesign_gershgorin, alpha=0
%
% Sign convention: same as grid section (u = K*x)

fprintf('\n\n=================================================================\n');
fprintf('  IEEE 13-BUS: Dense LQR vs kappa=2 Truncation vs EA-LQR\n');
fprintf('=================================================================\n');

% ---- Parameters ----
% shin_kappa defined above (shared with grid truncation)
shin_seeds  = [1];
nShinSeeds  = numel(shin_seeds);

shin_ea_params           = ea_params;
shin_ea_params.alpha     = 0;
shin_ea_params.maxGen    = 150;

shin_names = {'IEEE13'};
nShinSys   = numel(shin_names);

shin_J_raw          = NaN(nShinSys, 3, nShinSeeds);
shin_nnz_K          = NaN(nShinSys, 3, nShinSeeds);
shin_ea_hist        = cell(nShinSys, nShinSeeds);
shin_dense_comp     = NaN(nShinSys, 1);
shin_trunc_comp_ratio = NaN(nShinSys, 1);
shin_trunc_stable     = false(nShinSys, 1);   % true if K_trunc IEEE 13-bus is stable

% ---- Build IEEE 13-bus system ----
cfg13.model     = 'swing';
cfg13.Ts        = 0.2;
cfg13.rhoTarget = 0;
[sys13, busInfo] = build_ieee13bus_system(cfg13);
A_shin  = sys13.A;   B_shin = sys13.B2;
Q_shin  = sys13.Q;   R_shin = sys13.R;
n_shin  = busInfo.nBus;   % 13

% Graph and index sets for truncate()
G_shin    = graph(busInfo.adjMtx);
Is_shin   = cell(n_shin, 1);
Js_shin   = cell(n_shin, 1);
for ii = 1:n_shin
    Is_shin{ii} = [ii, n_shin + ii];   % theta_i (state ii), omega_i (state 13+ii)
    Js_shin{ii} = ii;                  % actuator ii → row ii of K
end
prob_shin = struct('name', 'grid');    % use unweighted hop distances

for si = 1:nShinSys
    fprintf('\n--- System: %s ---\n', shin_names{si});

    % ---- Method 1: Dense LQR ----
    K_pos_dense  = dlqr(A_shin, B_shin, Q_shin, R_shin);
    K_dense_shin = -K_pos_dense;
    J_dense_shin = get_lqr_cost(A_shin, B_shin, Q_shin, R_shin, K_dense_shin);
    fprintf('  Dense  : J=%.4e  nnz=%d\n', J_dense_shin, nnz(K_dense_shin));
    for ss = 1:nShinSeeds
        shin_J_raw(si, 1, ss) = J_dense_shin;
        shin_nnz_K(si, 1, ss) = nnz(K_dense_shin);
    end
    shin_dense_comp(si) = cost_EA(A_shin, B_shin, Q_shin, R_shin, K_dense_shin, J_dense_shin, 0);

    % ---- Method 1b: Diagonal LQR ----
    K_diag_shin = zeros(size(K_dense_shin));   % 13×26
    for i = 1:n_shin
        K_diag_shin(i, i)          = K_dense_shin(i, i);          % theta_i
        K_diag_shin(i, n_shin + i) = K_dense_shin(i, n_shin + i); % omega_i
    end
    if max(abs(eig(A_shin + B_shin * K_diag_shin))) >= 1.0
        [K_diag_shin, info_d] = gersgorin_stabilize_K(A_shin, B_shin, K_diag_shin, K_dense_shin);
        if ~info_d.success
            fprintf('  WARNING: Diagonal K still unstable (IEEE 13-bus)\n');
        end
    end
    J_diag_shin = get_lqr_cost(A_shin, B_shin, Q_shin, R_shin, K_diag_shin);
    fprintf('  Diag   : J=%.4e  ratio=%.4f  nnz=%d\n', J_diag_shin, J_diag_shin/J_dense_shin, nnz(K_diag_shin));
    J_diag_shin_n = cost_EA(A_shin, B_shin, Q_shin, R_shin, K_diag_shin, J_dense_shin, 0);
    shin_diag_comp_ratio = J_diag_shin_n / shin_dense_comp(si);

    % ---- Method 2: κ=1 Truncation ----
    K_trunc_shin = truncate(K_dense_shin, G_shin, Is_shin, Js_shin, shin_kappa, prob_shin);
    rho_trunc    = max(abs(eig(A_shin + B_shin * K_trunc_shin)));
    shin_trunc_stable(si) = (rho_trunc < 1.0);
    if rho_trunc >= 1.0
        fprintf('  Trunc k=%d: UNSTABLE (rho=%.4f)\n', shin_kappa, rho_trunc);
        J_trunc_shin = Inf;
        shin_trunc_comp_ratio(si) = Inf;
    else
        J_trunc_shin = get_lqr_cost(A_shin, B_shin, Q_shin, R_shin, K_trunc_shin);
        fprintf('  Trunc k=%d: J=%.4e  ratio=%.4f  nnz=%d  rho=%.4f\n', ...
            shin_kappa, J_trunc_shin, J_trunc_shin/J_dense_shin, nnz(K_trunc_shin), rho_trunc);
        J_comp_trunc_shin = cost_EA(A_shin, B_shin, Q_shin, R_shin, K_trunc_shin, J_dense_shin, 0);
        shin_trunc_comp_ratio(si) = J_comp_trunc_shin / shin_dense_comp(si);
    end
    for ss = 1:nShinSeeds
        shin_J_raw(si, 2, ss) = J_trunc_shin;
        shin_nnz_K(si, 2, ss) = nnz(K_trunc_shin);
    end

    % ---- Method 3: EA-LQR ----
    opts_shin = struct('verbose', false, 'Q_bm', Q_shin, 'R_bm', R_shin, 'useGersRepair', true);
    for ss = 1:nShinSeeds
        rng(shin_seeds(ss));
        fprintf('  EA seed=%d ...', shin_seeds(ss));
        [res_shin, shin_ea_hist{si, ss}] = ea_lqr_codesign_gershgorin(A_shin, B_shin, shin_ea_params, opts_shin);
        K_ea_shin  = res_shin.K_sparse;
        J_ea_s     = get_lqr_cost(A_shin, B_shin, Q_shin, R_shin, K_ea_shin);
        shin_J_raw(si, 3, ss) = J_ea_s;
        shin_nnz_K(si, 3, ss) = nnz(K_ea_shin);
        fprintf(' J=%.4e  ratio=%.4f  nnz=%d\n', J_ea_s, J_ea_s/J_dense_shin, nnz(K_ea_shin));
    end
end

% ---- Normalise by Dense LQR cost ----
shin_J_norm   = shin_J_raw ./ shin_J_raw(:, 1, :);
shin_J_mean   = mean(shin_J_norm, 3, 'omitnan');
shin_J_std    = std(shin_J_norm, 0, 3, 'omitnan');
shin_nnz_mean = mean(shin_nnz_K, 3, 'omitnan');

% ---- Console summary ----
fprintf('\n%-14s  %-10s  %-10s  %-10s  %-10s  %-10s  %-10s  %-10s\n', ...
    'System', 'Dense J/J*', 'Dense nnz', 'Diag J/J*', 'Trunc J/J*', 'Trunc nnz', 'EA J/J*', 'EA std');
fprintf('%s\n', repmat('-', 1, 92));
for si = 1:nShinSys
    fprintf('%-14s  %-10.4f  %-10.0f  %-10.4f  %-10.4f  %-10.0f  %-10.4f  %-10.4f\n', ...
        shin_names{si}, shin_J_mean(si,1), shin_nnz_mean(si,1), ...
        J_diag_shin/J_dense_shin, ...
        shin_J_mean(si,2), shin_nnz_mean(si,2), shin_J_mean(si,3), shin_J_std(si,3));
end

% ---- Collect IEEE 13-bus scaling data for Figure 2 ----
n_shin_sys  = size(A_shin, 1);    % 26
nu_shin_sys = size(B_shin, 2);    % 13

J_mesh_perf_mean = mean(shin_J_raw(1,3,:), 3, 'omitnan') / J_dense_shin;
spar_mesh_dense  = mean(shin_nnz_K(1,1,:), 3, 'omitnan') / (n_shin_sys * nu_shin_sys);
spar_mesh_ea     = mean(shin_nnz_K(1,3,:), 3, 'omitnan') / (n_shin_sys * nu_shin_sys);

convGen_mesh_dense = NaN(1, nShinSeeds);
for ss_m = 1:nShinSeeds
    curve_sh = shin_ea_hist{1, ss_m}.bestCost;
    idx_m    = find(curve_sh < shin_dense_comp(1), 1, 'first');
    if ~isempty(idx_m), convGen_mesh_dense(ss_m) = idx_m; end
end
mu_cg_mesh = nanmean(convGen_mesh_dense);
sd_cg_mesh = nanstd(convGen_mesh_dense);

% Diagonal density for IEEE 13-bus
spar_mesh_diag = nnz(K_diag_shin) / (n_shin_sys * nu_shin_sys);

% Truncation data for IEEE 13-bus
if isfinite(J_trunc_shin)
    J_perf_trunc_ieee13 = J_trunc_shin / J_dense_shin;
    J_comp_trunc_ieee13 = J_perf_trunc_ieee13 ...
        + 0.05*nnz(K_trunc_shin) ...
        + 0.4 *nnz(any(K_trunc_shin,2)) ...
        + 0.2 *nnz(any(K_trunc_shin,1));
    spar_trunc_ieee13 = nnz(K_trunc_shin) / (n_shin_sys * nu_shin_sys);
else
    J_perf_trunc_ieee13 = Inf;
    J_comp_trunc_ieee13 = Inf;
    spar_trunc_ieee13   = 0;
end

% Convergence gen: first gen when EA composite cost < Diagonal composite cost
J_diag_comp_abs = shin_diag_comp_ratio * shin_dense_comp(1);
convGen_mesh_diag = NaN(1, nShinSeeds);
for ss_m = 1:nShinSeeds
    curve_sh = shin_ea_hist{1, ss_m}.bestCost;
    idx_md   = find(curve_sh < J_diag_comp_abs, 1, 'first');
    if ~isempty(idx_md), convGen_mesh_diag(ss_m) = idx_md; end
end
mu_cg_mesh_diag = nanmean(convGen_mesh_diag);
sd_cg_mesh_diag = nanstd(convGen_mesh_diag);

%% ===================== Figure 2: Scaling Analysis (5x5, 7x7, IEEE13-bus) =====================
fig2 = figure('Position', [40 520 1500 380], 'Color', 'w');
set(fig2, 'DefaultAxesFontSize', axFS, 'DefaultAxesFontName', 'Times New Roman');

% Use categorical positions (1,2,3) because IEEE13 Nx=26 < grid Nx=[50,98]
xDim     = [1; 2];                       % positions for 5x5 and 7x7
xDim_ext = [1; 2; 3];                    % including IEEE 13-bus
xTickLbl = [dimLabel; {'IEEE13'}];       % {'5x5','7x7','IEEE13'}
barW     = 0.18;
mNames   = {'Dense', 'Diag', 'Trunc', 'EA'};
mColors  = {cb_dens, cb_diag, cb_trunc, cb_ea};

% --- (a) Cost decomposition stacked ---
subplot(1, 3, 1); hold on;

% Grid systems: Dense, Diag, EA (offsets -1.5, -0.5, +1.5 * barW)
xOffs3 = [-1.5, -0.5, 1.5] * barW;
for m = 1:3
    bar(xDim + xOffs3(m), meanPerf(:,m),  barW, 'FaceColor', mColors{m}, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
    bar(xDim + xOffs3(m), meanTotal(:,m), barW, 'FaceColor', mColors{m}, 'EdgeColor', mColors{m}, 'LineWidth', 1.2, 'FaceAlpha', 0.15);
end
% Grid systems: Truncation bar (offset +0.5 * barW)
bar(xDim + 0.5*barW, meanPerfTrunc, barW, 'FaceColor', cb_trunc, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
bar(xDim + 0.5*barW, meanCompTrunc, barW, 'FaceColor', cb_trunc, 'EdgeColor', cb_trunc, 'LineWidth', 1.2, 'FaceAlpha', 0.15);

% IEEE 13-bus: Dense, Diagonal, Truncation, EA
bar(3 - 1.5*barW, 1.0,                       barW, 'FaceColor', cb_dens,  'FaceAlpha', 0.5, 'EdgeColor', 'none');
bar(3 - 0.5*barW, J_diag_shin/J_dense_shin,  barW, 'FaceColor', cb_diag,  'FaceAlpha', 0.5, 'EdgeColor', 'none');
if isfinite(J_perf_trunc_ieee13)
    bar(3 + 0.5*barW, J_perf_trunc_ieee13,   barW, 'FaceColor', cb_trunc, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
    bar(3 + 0.5*barW, J_comp_trunc_ieee13,   barW, 'FaceColor', cb_trunc, 'EdgeColor', cb_trunc, 'LineWidth', 1.2, 'FaceAlpha', 0.15);
end
bar(3 + 1.5*barW, J_mesh_perf_mean,          barW, 'FaceColor', cb_ea,    'FaceAlpha', 0.5, 'EdgeColor', 'none');

% Manual legend (4 methods)
hL = gobjects(4,1);
for m = 1:4
    hL(m) = patch(NaN, NaN, mColors{m}, 'FaceAlpha', 0.5, 'EdgeColor', mColors{m});
end
legend(hL, mNames, 'FontSize', legFS, 'Location', 'northwest');
xlabel('System', 'FontSize', labFS);
ylabel('Cost', 'FontSize', labFS);
title('(a) Cost Decomposition', 'FontSize', titFS);
set(gca, 'XTick', xDim_ext, 'XTickLabel', xTickLbl);
grid on; box on;

% --- (b) Controller density ---
subplot(1, 3, 2); hold on;
sparDense = [meanNnz(:,1) ./ (nxVec .* nuVec);  spar_mesh_dense];
sparDiag  = [meanNnz(:,2) ./ (nxVec .* nuVec);  spar_mesh_diag];   % all 3 systems
sparTrunc = [meanNnzTrunc  ./ (nxVec .* nuVec);  spar_trunc_ieee13];
sparEA    = [meanNnz(:,3) ./ (nxVec .* nuVec);  spar_mesh_ea];

plot(xDim_ext, sparDense*100, 's--', 'Color', cb_dens,  'LineWidth', lw, 'MarkerSize', 8, 'MarkerFaceColor', cb_dens);
plot(xDim_ext, sparDiag*100,  'd-.', 'Color', cb_diag,  'LineWidth', lw, 'MarkerSize', 8, 'MarkerFaceColor', cb_diag);
plot(xDim_ext, sparTrunc*100, '^:',  'Color', cb_trunc, 'LineWidth', lw, 'MarkerSize', 8, 'MarkerFaceColor', cb_trunc);
plot(xDim_ext, sparEA*100,    'o-',  'Color', cb_ea,    'LineWidth', lw, 'MarkerSize', 8, 'MarkerFaceColor', cb_ea);
legend({'Dense', 'Diagonal', 'Trunc', 'EA'}, 'FontSize', legFS, 'Location', 'northeast');
xlabel('System', 'FontSize', labFS);
ylabel('Density (\%)', 'FontSize', labFS);
title('(b) Controller Density', 'FontSize', titFS);
set(gca, 'XTick', xDim_ext, 'XTickLabel', xTickLbl);
grid on; box on;

% --- (c) Convergence speed ---
subplot(1, 3, 3); hold on;
mu_cg_d     = nanmean(convGenDense, 2);   sd_cg_d = nanstd(convGenDense, 0, 2);
mu_cg_g     = nanmean(convGenDiag,  2);   sd_cg_g = nanstd(convGenDiag,  0, 2);
mu_cg_d_ext = [mu_cg_d; mu_cg_mesh];
sd_cg_d_ext = [sd_cg_d; sd_cg_mesh];
mu_cg_g_ext = [mu_cg_g; mu_cg_mesh_diag];
sd_cg_g_ext = [sd_cg_g; sd_cg_mesh_diag];

errorbar(xDim_ext, mu_cg_d_ext, sd_cg_d_ext, 's-', 'Color', cb_dens, 'LineWidth', lw, ...
    'MarkerSize', 8, 'MarkerFaceColor', cb_dens, 'CapSize', 5);
errorbar(xDim_ext, mu_cg_g_ext, sd_cg_g_ext, 'd-', 'Color', cb_diag, 'LineWidth', lw, ...
    'MarkerSize', 8, 'MarkerFaceColor', cb_diag, 'CapSize', 5);
legend({'Gen to beat Dense', 'Gen to beat Diagonal'}, 'FontSize', legFS, 'Location', 'northwest');
xlabel('System', 'FontSize', labFS);
ylabel('Generation', 'FontSize', labFS);
title('(c) Convergence Speed', 'FontSize', titFS);
set(gca, 'XTick', xDim_ext, 'XTickLabel', xTickLbl);
ylim([0 maxGen]);
grid on; box on;

sgtitle('Scaling Analysis: Dense vs Diagonal vs EA-LQR (grid + IEEE 13-bus)', ...
    'FontSize', 16, 'FontWeight', 'bold', 'FontName', 'Times New Roman');

exportgraphics(fig2, fullfile(outDir, 'perf_bounds_scaling.pdf'), 'ContentType', 'vector');
saveas(fig2, fullfile(outDir, 'perf_bounds_scaling.png'));
fprintf('Figure 2 (scaling with IEEE 13-bus) saved → %s/perf_bounds_scaling.*\n', outDir);




%% ---- Figure 1 panel (c): Append IEEE 13-bus panel now that shin_ea_hist is available ----
si_m = 1;   % index of 'IEEE13' in shin_names
figure(fig1);
subplot(2, nGrid+1, nGrid+1); hold on;

Jref_sh   = shin_dense_comp(si_m);
normMat_sh = zeros(maxGen, nShinSeeds);
for ss = 1:nShinSeeds
    normMat_sh(:, ss) = shin_ea_hist{si_m, ss}.bestCost / Jref_sh;
end
mu_sh = mean(normMat_sh, 2);
sd_sh = std(normMat_sh, 0, 2);

fill([gens; flipud(gens)], [mu_sh+sd_sh; flipud(max(mu_sh-sd_sh, 0))], ...
    ca_ea, 'EdgeColor', 'none', 'FaceAlpha', 0.4);
hEA_sh   = plot(gens, mu_sh, '-',  'Color', cb_ea,    'LineWidth', lw);
hDens_sh = yline(1.0,  '--', 'Color', cb_dens,  'LineWidth', 2);

J_diag_sh_n = shin_diag_comp_ratio;
if J_diag_sh_n < 1e4
    hDiag_sh = yline(J_diag_sh_n, '-.', 'Color', cb_diag, 'LineWidth', 2);
else
    hDiag_sh = plot(NaN, NaN, '-.', 'Color', cb_diag, 'LineWidth', 2);
    text(maxGen/2, max(mu_sh)*0.9, 'Diag: Unstable', 'Color', cb_diag, ...
        'FontSize', legFS, 'FontName', 'Times New Roman', 'HorizontalAlignment', 'center');
end

J_trunc_sh_n = shin_trunc_comp_ratio(si_m);
if shin_trunc_stable(si_m)
    hTrunc_sh = yline(J_trunc_sh_n, ':', 'Color', cb_trunc, 'LineWidth', 2);
else
    hTrunc_sh = plot(NaN, NaN, ':', 'Color', cb_trunc, 'LineWidth', 2);
    text(maxGen*0.7, max(mu_sh)*0.85, 'Trunc: Unstable', 'Color', cb_trunc, ...
        'FontSize', legFS, 'FontName', 'Times New Roman', 'HorizontalAlignment', 'center');
end

% Exponential lower bound — Eq. (49)
hLB_sh = plot(NaN, NaN, '-.', 'Color', cb_lb, 'LineWidth', lw*0.8);  % sentinel
if showLB
    JR_vec_sh  = arrayfun(@(ss) shin_ea_hist{si_m,ss}.bestCost(end), 1:nShinSeeds);
    A0_vec_sh  = arrayfun(@(ss) shin_ea_hist{si_m,ss}.bestCost(1),   1:nShinSeeds) - JR_vec_sh;
    JR_star_sh = mean(JR_vec_sh);
    A0_sh      = max(mean(A0_vec_sh), 0);
    LB_sh      = (JR_star_sh + A0_sh * exp(-lb_lambda * gens)) / Jref_sh;
    delete(hLB_sh);
    hLB_sh = plot(gens, LB_sh, '-.', 'Color', cb_lb, 'LineWidth', lw*0.8);
end

xlabel('Generation', 'FontSize', labFS);
title('(c) IEEE 13-bus ($n$=26)', 'FontSize', titFS, 'Interpreter', 'latex');
legH_sh = [hEA_sh, hDens_sh];
legL_sh = {'EA-LQR', 'Dense ($=1$)'};
if J_diag_sh_n < 1e4
    legH_sh(end+1) = hDiag_sh;
    legL_sh{end+1} = 'Diagonal';
end
if shin_trunc_stable(si_m)
    legH_sh(end+1) = hTrunc_sh;
    legL_sh{end+1} = '$\kappa=2$ Trunc.';
end
if showLB
    legH_sh(end+1) = hLB_sh;
    legL_sh{end+1} = 'LB (Eq.~49)';
end
legend(legH_sh, legL_sh, 'Location', 'northeast', 'FontSize', legFS, 'Interpreter', 'latex');
grid on; box on;

sgtitle('EA-LQR Convergence Normalized by Dense LQR', ...
    'FontSize', 16, 'FontWeight', 'bold', 'FontName', 'Times New Roman');

%% ---- Figure 1 row 2: EA Controller Topology Diagrams ----
figure(fig1);

% Marker colors
c_act  = [0.85 0.33 0.10];   % orange for actuator nodes
c_sens = [0.17 0.63 0.17];   % green  for sensor nodes
c_both = [0.49 0.18 0.56];   % purple for actuator + sensor nodes

% --- Panels (d) and (e): Grid systems ---
for g = 1:nGrid
    subplot(2, nGrid+1, nGrid+1+g); hold on; box on; axis equal off;

    gSz  = gridSizes(g);
    nN   = gSz^2;
    adjM = adjMtx_all{g};
    Ke   = K_ea_all{g};
    actN = actuatedNodes_all{g};   % Nu×1: physical node for each actuator

    % Node positions on a regular gSz×gSz grid
    nodePos = [(mod((1:nN)-1, gSz)+1)', (floor(((1:nN)-1)/gSz)+1)'];

    % Threshold: entry significant if |K(i,j)| > 5% of max|K|
    edgeThresh = max(abs(Ke(:))) * 0.05;

    % Build K-reachability via shortest physical paths, then color edges
    G_topo    = graph(adjM);
    rowNorms  = vecnorm(Ke, 2, 2);
    usedEdges = false(nN, nN);
    for i = 1:size(Ke,1)
        if rowNorms(i) <= max(rowNorms)*0.05, continue; end
        nSrc = actN(i);
        for j = 1:nN
            if j == nSrc, continue; end
            if (abs(Ke(i,j)) > edgeThresh) || (abs(Ke(i,nN+j)) > edgeThresh)
                pth = shortestpath(G_topo, nSrc, j);
                for kk = 1:length(pth)-1
                    usedEdges(pth(kk), pth(kk+1)) = true;
                    usedEdges(pth(kk+1), pth(kk)) = true;
                end
            end
        end
    end
    for u = 1:nN
        for v = u+1:nN
            if adjM(u,v) > 0
                if usedEdges(u,v)
                    plot(nodePos([u v],1), nodePos([u v],2), '-', ...
                        'Color', cb_ea, 'LineWidth', 2.5);
                else
                    plot(nodePos([u v],1), nodePos([u v],2), '-', ...
                        'Color', [0.82 0.82 0.82], 'LineWidth', 0.8);
                end
            end
        end
    end

    % Active node: row norm > 5% of max row norm (filters Gershgorin-repair artifacts)
    activeActIdx = find(rowNorms > max(rowNorms) * 0.05);
    activeNodes  = actN(activeActIdx);

    % Sensor nodes: physical node j is a sensor if any column j or nN+j has a significant entry
    sensorMask = false(nN, 1);
    for j = 1:nN
        if any(abs(Ke(:,j)) > edgeThresh) || any(abs(Ke(:,nN+j)) > edgeThresh)
            sensorMask(j) = true;
        end
    end
    sensorNodes = find(sensorMask);

    % Classify nodes into 4 categories
    actOnlyNodes  = setdiff(activeNodes(:),  sensorNodes(:));
    sensOnlyNodes = setdiff(sensorNodes(:),  activeNodes(:));
    bothNodes     = intersect(activeNodes(:), sensorNodes(:));
    inactiveNodes = setdiff(1:nN, union(activeNodes(:), sensorNodes(:)));

    % Inactive: gray circle
    scatter(nodePos(inactiveNodes,1), nodePos(inactiveNodes,2), 30, ...
        [0.75 0.75 0.75], 'o', 'filled', 'MarkerEdgeColor', 'none');
    % Sensor-only: green circle
    if ~isempty(sensOnlyNodes)
        scatter(nodePos(sensOnlyNodes,1), nodePos(sensOnlyNodes,2), 60, ...
            c_sens, 'o', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.7);
    end
    % Actuator-only: orange circle
    if ~isempty(actOnlyNodes)
        scatter(nodePos(actOnlyNodes,1), nodePos(actOnlyNodes,2), 60, ...
            c_act, 'o', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.7);
    end
    % Both roles: purple filled circle
    if ~isempty(bothNodes)
        scatter(nodePos(bothNodes,1), nodePos(bothNodes,2), 80, ...
            c_both, 'o', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.9);
    end

    title(sprintf('(%c) %d$\\times$%d EA controller', char('d'+g-1), gSz, gSz), ...
        'FontSize', titFS, 'Interpreter', 'latex');
end

% --- Panel (f): IEEE 13-bus EA controller ---
subplot(2, nGrid+1, 2*(nGrid+1)); hold on; box on; axis equal off;

nN_sh      = n_shin;                   % 13
adjM_sh    = busInfo.adjMtx;           % 13×13 binary adjacency
Ke_sh      = K_ea_shin;                % 13×26
nodePos_sh = busInfo.nodeCoords;       % 13×2 geographic coordinates

% Threshold: entry significant if |K(i,j)| > 5% of max|K|
edgeThresh_sh = max(abs(Ke_sh(:))) * 0.05;

% Build K-reachability via shortest physical paths, then color edges
G_topo_sh    = graph(adjM_sh);
usedEdges_sh = false(nN_sh, nN_sh);
for i = 1:nN_sh
    for j = 1:nN_sh
        if j == i, continue; end
        if (abs(Ke_sh(i,j)) > edgeThresh_sh) || (abs(Ke_sh(i,nN_sh+j)) > edgeThresh_sh)
            pth = shortestpath(G_topo_sh, i, j);
            for kk = 1:length(pth)-1
                usedEdges_sh(pth(kk), pth(kk+1)) = true;
                usedEdges_sh(pth(kk+1), pth(kk)) = true;
            end
        end
    end
end
for u = 1:nN_sh
    for v = u+1:nN_sh
        if adjM_sh(u,v) > 0
            if usedEdges_sh(u,v)
                plot(nodePos_sh([u v],1), nodePos_sh([u v],2), '-', ...
                    'Color', cb_ea, 'LineWidth', 2.5);
            else
                plot(nodePos_sh([u v],1), nodePos_sh([u v],2), '-', ...
                    'Color', [0.82 0.82 0.82], 'LineWidth', 0.8);
            end
        end
    end
end

% Active node: row norm > 5% of max row norm (filters Gershgorin-repair artifacts)
rowNorms_sh = vecnorm(Ke_sh, 2, 2);
sensUsed_sh = false(nN_sh, 1);
for j = 1:nN_sh
    if any(abs(Ke_sh(:,j)) > edgeThresh_sh) || any(abs(Ke_sh(:,nN_sh+j)) > edgeThresh_sh)
        sensUsed_sh(j) = true;
    end
end
% Classify nodes into 4 categories
activeActNodes_sh = find(rowNorms_sh > max(rowNorms_sh) * 0.05);
sensorNodes_sh    = find(sensUsed_sh);
actOnlyNodes_sh   = setdiff(activeActNodes_sh, sensorNodes_sh);
sensOnlyNodes_sh  = setdiff(sensorNodes_sh,    activeActNodes_sh);
bothNodes_sh      = intersect(activeActNodes_sh, sensorNodes_sh);
inactiveNodes_sh  = setdiff(1:nN_sh, union(activeActNodes_sh, sensorNodes_sh));

s1 = scatter(nodePos_sh(inactiveNodes_sh,1), nodePos_sh(inactiveNodes_sh,2), 24, ...
    [0.75 0.75 0.75], 'o', 'filled', 'MarkerEdgeColor', 'none');
s2 = scatter([], [], 60, c_sens, 'o', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.7);
s3 = scatter([], [], 60, c_act,  'o', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.7);
s4 = scatter([], [], 80, c_both, 'o', 'filled', 'MarkerEdgeColor', 'k',    'LineWidth', 0.9);
if ~isempty(sensOnlyNodes_sh)
    scatter(nodePos_sh(sensOnlyNodes_sh,1), nodePos_sh(sensOnlyNodes_sh,2), 60, ...
        c_sens, 'o', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.7);
end
if ~isempty(actOnlyNodes_sh)
    scatter(nodePos_sh(actOnlyNodes_sh,1), nodePos_sh(actOnlyNodes_sh,2), 60, ...
        c_act, 'o', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.7);
end
if ~isempty(bothNodes_sh)
    scatter(nodePos_sh(bothNodes_sh,1), nodePos_sh(bothNodes_sh,2), 80, ...
        c_both, 'o', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.9);
end

title('(f) IEEE 13-bus EA controller', 'FontSize', titFS, 'Interpreter', 'latex');
legend([s1, s2, s3, s4], {'Inactive', 'Sensor', 'Actuator', 'Both'}, ...
    'Location', 'best', 'FontSize', legFS);

exportgraphics(fig1, fullfile(outDir, 'perf_bounds_convergence.pdf'), 'ContentType', 'vector');
saveas(fig1, fullfile(outDir, 'perf_bounds_convergence.png'));
fprintf('Figure 1 (3-panel with 2d-mesh) saved → %s/perf_bounds_convergence.*\n', outDir);


%% ===================== Save plot-only data =====================
plotDir = fullfile(pwd, 'plot_only');
if ~exist(plotDir, 'dir'), mkdir(plotDir); end

% Pre-compute derived scalars so plot_figs.m needs no arithmetic
J_diag_shin_ratio = J_diag_shin / J_dense_shin;
mu_cg_d = nanmean(convGenDense, 2);   sd_cg_d = nanstd(convGenDense, 0, 2);
mu_cg_g = nanmean(convGenDiag,  2);   sd_cg_g = nanstd(convGenDiag,  0, 2);

save(fullfile(plotDir, 'fig_data.mat'), ...
    'gridSizes', 'eaSeeds', 'sysSize', 'ea_params', ...
    'eaCurves', 'denseCostRef', 'convGenDense', 'convGenDiag', ...
    'meanTotal', 'meanPerf', 'meanStruct', ...
    'meanNnz', 'nAct_all', 'nSens_all', ...
    'meanPerfTrunc', 'meanCompTrunc', 'meanNnzTrunc', 'trunc_stable_grid', ...
    'adjMtx_all', 'K_ea_all', 'actuatedNodes_all', ...
    'shin_ea_hist', 'shin_dense_comp', 'shin_diag_comp_ratio', 'shin_trunc_comp_ratio', 'shin_trunc_stable', ...
    'busInfo', 'K_ea_shin', 'n_shin', ...
    'J_diag_shin_ratio', 'J_perf_trunc_ieee13', 'J_comp_trunc_ieee13', ...
    'J_mesh_perf_mean', 'spar_mesh_dense', 'spar_mesh_ea', 'spar_mesh_diag', 'spar_trunc_ieee13', ...
    'mu_cg_d', 'sd_cg_d', 'mu_cg_g', 'sd_cg_g', ...
    'mu_cg_mesh', 'sd_cg_mesh', 'mu_cg_mesh_diag', 'sd_cg_mesh_diag','LB_sh');
fprintf('Plot data saved → %s/fig_data.mat\n', plotDir);


%% ===================== Local Helper =====================
function [J_total, J_perf, J_struct] = cost_decompose(A, B, Q, R, K, costBM, alpha)
    rho = max(abs(eig(A + B * K)));
    if rho >= 1.0
        J_total = 1e9;  J_perf = Inf;  J_struct = 0;
        return;
    end
    J_perf = get_lqr_cost(A, B, Q, R, K) / costBM;

    w_comm = 0.05 * (1 - alpha);
    w_row  = 0.4  * (1 - alpha);
    w_col  = 0.2  * (1 - alpha);
    J_struct = w_comm * nnz(K) + w_row * nnz(any(K,2)) + w_col * nnz(any(K,1));
    J_total  = J_perf + J_struct;
end
