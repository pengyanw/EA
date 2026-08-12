%% plot_figs.m
% Reproduce Fig 1 (convergence + topology) and Fig 2 (unstable: Baseline vs Gershgorin)
% from saved data — no algorithm required.
%
% Prerequisites (run once each to generate the .mat files):
%   analysis_perf_bounds.m    → fig_data.mat
%   run_multi_seed_unstable.m → fig_data_unstable.mat
%
% Usage:
%   cd plot_only
%   plot_figs
%
% This script only draws; every curve is read from the .mat files. LB_grid / LB_sh
% now hold the Theorem C' certified LOWER BOUND on min_a J_EA (computed once in
% analysis_perf_bounds.m via EA functions/gap_predictor.m), not the superseded
% Theorem 1 convergence prediction of Phi_predictor.m -- so there is no second
% copy of the formula here to drift out of sync.

clear; clc; close all;
load('fig_data.mat');
U = load('fig_data_unstable.mat');   % isolated struct — avoids overwriting numSeeds, eaSeeds, etc.

%% ===================== Configurable display options =====================
showLB = true;   % Show the Thm. C' certified optimality gap (Section IV)

% shin_kappa and LB_grid were added to fig_data.mat later than the other fields;
% fall back gracefully so an older .mat still plots.
if ~exist('shin_kappa', 'var'), shin_kappa = 1; end

%% ===================== Derive simple quantities =====================
nxVec    = sysSize(:,1);
nuVec    = sysSize(:,2);
nGrid    = numel(gridSizes);
numSeeds = numel(eaSeeds);
maxGen   = ea_params.maxGen;
gens     = (1:maxGen)';
dimLabel = arrayfun(@(g) sprintf('%dx%d', g, g), gridSizes, 'UniformOutput', false)';
nShinSeeds = size(shin_ea_hist, 2);

%% ===================== Colors & fonts =====================
axFS = 18;  labFS = 19;  titFS = 19;  legFS = 15;  lw = 3;
cb_ea     = [0.00 0.45 0.74];
cb_dens   = [0.50 0.50 0.50];
cb_diag   = [0.85 0.33 0.10];
cb_trunc  = [0.47 0.67 0.19];
cb_gd     = [0.49 0.18 0.56];   % purple — greedy baseline (B2 comparison)
cb_lb     = [0.93 0.69 0.13];   % amber — Thm. C' certified lower envelope
ca_lb     = [0.99 0.92 0.72];   % amber fill — the certified optimality gap
cb_purple = [0.49 0.18 0.56];   % purple — act + sensor node
ca_ea     = [0.75 0.85 1.00];
barW     = 0.18;

%% ===================== Figure 1: Normalized Convergence =====================
fig1 = figure('Position', [40 60 2200 820], 'Color', 'w');
set(fig1, 'DefaultAxesFontSize', axFS, 'DefaultAxesFontName', 'Times New Roman');

for g = 1:nGrid
    subplot(2, nGrid+1, g); hold on;

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

    % Greedy baseline on the same budget axis (its evaluations / N_p).
    hGD = gobjects(0);
    if exist('gdCurves', 'var') && ~isempty(gdCurves{g, 1})
        gdMat = zeros(maxGen, numSeeds);
        for s = 1:numSeeds
            gdMat(:, s) = gdCurves{g, s} / denseCostRef(g, s);
        end
        hGD = plot(gens, mean(gdMat, 2), '-', 'Color', cb_gd, 'LineWidth', lw*0.8);
    end

    hDens = yline(1.0, '--', 'Color', cb_dens, 'LineWidth', 2);

    J_diag_n = meanTotal(g,2) / Jref;
    if J_diag_n < 1e4
        hDiag = yline(J_diag_n, '-.', 'Color', cb_diag, 'LineWidth', 2);
    else
        hDiag = plot(NaN, NaN, '-.', 'Color', cb_diag, 'LineWidth', 2);
        text(maxGen/2, max(mu)*0.9, 'Diag: Unstable', 'Color', cb_diag, ...
            'FontSize', legFS, 'FontName', 'Times New Roman', 'HorizontalAlignment', 'center');
    end

    J_trunc_n = meanCompTrunc(g) / Jref;
    if all(trunc_stable_grid(g, :))
        hTrunc = yline(J_trunc_n, ':', 'Color', cb_trunc, 'LineWidth', 2);
    else
        hTrunc = plot(NaN, NaN, ':', 'Color', cb_trunc, 'LineWidth', 2);
        text(maxGen*0.7, max(mu)*0.85, 'Trunc: Unstable', 'Color', cb_trunc, ...
            'FontSize', legFS, 'FontName', 'Times New Roman', 'HorizontalAlignment', 'center');
    end

    % Theorem C' certified lower bound on min_a J_EA (precomputed, see header).
    % Shade between it and the EA curve: that band is the certified optimality gap.
    hLB = plot(NaN, NaN, '-.', 'Color', cb_lb, 'LineWidth', lw*0.8);
    if showLB && exist('LB_grid', 'var') && any(isfinite(LB_grid(:, g)))
        delete(hLB);
        fill([gens; flipud(gens)], [mu; flipud(LB_grid(:, g))], ca_lb, ...
            'EdgeColor', 'none', 'FaceAlpha', 0.55, 'HandleVisibility', 'off');
        hLB = plot(gens, LB_grid(:, g), '-.', 'Color', cb_lb, 'LineWidth', lw*0.8);
        uistack(hEA, 'top');
    end

    % Mark convergence generation (mean)
    cg = round(nanmean(convGenDense(g,:)));
    if ~isnan(cg) && cg <= maxGen
        plot(cg, mu(cg), 'v', 'MarkerSize', 8, 'MarkerFaceColor', cb_ea, ...
            'MarkerEdgeColor', 'k', 'LineWidth', 0.8);
    end

    % Common budget axis: the EA spends N_p evaluations per generation, and the
    % greedy incumbent is sampled at N_p*g, so both curves share this scale.
    xlabel('Generation ($=$ evals$/N_p$)', 'FontSize', labFS, 'Interpreter', 'latex');
    if g == 1
        ylabel('$J \;/\; J_{\mathrm{dense}}$', 'FontSize', labFS, 'Interpreter', 'latex');
    end
    title(sprintf('(%c) %d$\\times$%d ($n$=%d, $m$=%d)', ...
        char('a'+g-1), gridSizes(g), gridSizes(g), nxVec(g), nuVec(g)), ...
        'FontSize', titFS, 'Interpreter', 'latex');
    grid on; box on;

    if g == 1 || g == nGrid
        legHandles = hEA;  legLabels = {'EA-LQR'};
        if ~isempty(hGD), legHandles(end+1) = hGD; legLabels{end+1} = 'Greedy'; end
        legHandles = [legHandles, hDens, hDiag];
        legLabels  = [legLabels, {'Dense ($=1$)', 'Diagonal'}];
        if all(trunc_stable_grid(g, :))
            legHandles(end+1) = hTrunc;
            legLabels{end+1}  = sprintf('$\\kappa=%d$ Trunc.', shin_kappa);
        end
        if showLB
            legHandles(end+1) = hLB;
            legLabels{end+1}  = 'Thm.~C bound on $\min_{\mathbf{a}} J_{\mathrm{EA}}$';
        end
        legend(legHandles, legLabels, 'Location', 'northeast', ...
            'FontSize', legFS, 'Interpreter', 'latex');
    end
end

sgtitle('EA-LQR Convergence Normalized by Dense LQR', ...
    'FontSize', 16, 'FontWeight', 'bold', 'FontName', 'Times New Roman');

%% ---- Figure 1 panel (c): IEEE 13-bus ----
si_m = 1;
figure(fig1);
subplot(2, nGrid+1, nGrid+1); hold on;

Jref_sh    = shin_dense_comp(si_m);
normMat_sh = zeros(maxGen, nShinSeeds);
for ss = 1:nShinSeeds
    normMat_sh(:, ss) = shin_ea_hist{si_m, ss}.bestCost / Jref_sh;
end
mu_sh = mean(normMat_sh, 2);
sd_sh = std(normMat_sh, 0, 2);

fill([gens; flipud(gens)], [mu_sh+sd_sh; flipud(max(mu_sh-sd_sh, 0))], ...
    ca_ea, 'EdgeColor', 'none', 'FaceAlpha', 0.4);
hEA_sh   = plot(gens, mu_sh, '-',  'Color', cb_ea,   'LineWidth', lw);
hGD_sh   = gobjects(0);
if exist('gdCurve_sh', 'var') && ~isempty(gdCurve_sh)
    hGD_sh = plot(gens, gdCurve_sh / Jref_sh, '-', 'Color', cb_gd, 'LineWidth', lw*0.8);
end
hDens_sh = yline(1.0, '--', 'Color', cb_dens, 'LineWidth', 2);

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
hLB_sh = plot(NaN, NaN, '-.', 'Color', cb_lb, 'LineWidth', lw*0.8);
if showLB && exist('LB_sh', 'var') && any(isfinite(LB_sh))
    delete(hLB_sh);
    fill([gens; flipud(gens)], [mu_sh; flipud(LB_sh)], ca_lb, ...
        'EdgeColor', 'none', 'FaceAlpha', 0.55, 'HandleVisibility', 'off');
    hLB_sh = plot(gens, LB_sh, '-.', 'Color', cb_lb, 'LineWidth', lw*0.8);
    uistack(hEA_sh, 'top');
end
xlabel('Generation ($=$ evals$/N_p$)', 'FontSize', labFS, 'Interpreter', 'latex');
title('(c) IEEE 13-bus ($n$=26)', 'FontSize', titFS, 'Interpreter', 'latex');
legH_sh = hEA_sh;  legL_sh = {'EA-LQR'};
if ~isempty(hGD_sh), legH_sh(end+1) = hGD_sh; legL_sh{end+1} = 'Greedy'; end
legH_sh(end+1) = hDens_sh;  legL_sh{end+1} = 'Dense ($=1$)';
if J_diag_sh_n < 1e4
    legH_sh(end+1) = hDiag_sh;
    legL_sh{end+1} = 'Diagonal';
end
if shin_trunc_stable(si_m)
    legH_sh(end+1) = hTrunc_sh;
    legL_sh{end+1} = sprintf('$\\kappa=%d$ Trunc.', shin_kappa);
end
if showLB
    legH_sh(end+1) = hLB_sh;
    legL_sh{end+1} = 'Thm.~C bound on $\min_{\mathbf{a}} J_{\mathrm{EA}}$';
end
legend(legH_sh, legL_sh, 'Location', 'northeast', 'FontSize', legFS, 'Interpreter', 'latex');
grid on; box on;

sgtitle('EA-LQR Convergence Normalized by Dense LQR', ...
    'FontSize', 16, 'FontWeight', 'bold', 'FontName', 'Times New Roman');

%% ---- Figure 1 row 2: EA Controller Topology Diagrams ----
figure(fig1);

% --- Panels (d) and (e): Grid systems ---
for g = 1:nGrid
    subplot(2, nGrid+1, nGrid+1+g); hold on; box on; axis equal off;

    gSz  = gridSizes(g);
    nN   = gSz^2;
    adjM = adjMtx_all{g};
    Ke   = K_ea_all{g};
    actN = actuatedNodes_all{g};

    nodePos = [(mod((1:nN)-1, gSz)+1)', (floor(((1:nN)-1)/gSz)+1)'];
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

    activeActIdx = find(rowNorms > max(rowNorms) * 0.05);
    actSet       = actN(activeActIdx);            % active actuator nodes

    sensFlags = false(nN, 1);
    for jj = 1:nN
        if any(abs(Ke(:,jj)) > edgeThresh) || any(abs(Ke(:,nN+jj)) > edgeThresh)
            sensFlags(jj) = true;
        end
    end
    sensSet = find(sensFlags);

    bothSet     = intersect(actSet(:), sensSet(:));
    actOnlySet  = setdiff(actSet(:),  bothSet);
    sensOnlySet = setdiff(sensSet(:), bothSet);
    inactiveSet = setdiff((1:nN)', union(actSet(:), sensSet(:)));

    scatter(nodePos(inactiveSet ,1), nodePos(inactiveSet ,2), 36, [0.75 0.75 0.75], 'filled', 'MarkerEdgeColor', 'none');
    scatter(nodePos(sensOnlySet ,1), nodePos(sensOnlySet ,2), 56, cb_trunc,  'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.7);
    scatter(nodePos(actOnlySet  ,1), nodePos(actOnlySet  ,2), 56, cb_diag,   'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.7);
    scatter(nodePos(bothSet     ,1), nodePos(bothSet     ,2), 70, cb_purple, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.9);

    title(sprintf('(%c) %d$\\times$%d EA controller', char('d'+g-1), gSz, gSz), ...
        'FontSize', titFS, 'Interpreter', 'latex');
end

% --- Panel (f): IEEE 13-bus EA controller ---
subplot(2, nGrid+1, 2*(nGrid+1)); hold on; box on; axis equal off;

nN_sh      = n_shin;
adjM_sh    = busInfo.adjMtx;
Ke_sh      = K_ea_shin;
nodePos_sh = busInfo.nodeCoords;

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

rowNorms_sh = vecnorm(Ke_sh, 2, 2);
sensUsed_sh = false(nN_sh, 1);
for j = 1:nN_sh
    if any(abs(Ke_sh(:,j)) > edgeThresh_sh) || any(abs(Ke_sh(:,nN_sh+j)) > edgeThresh_sh)
        sensUsed_sh(j) = true;
    end
end
actSet_sh      = find(rowNorms_sh > max(rowNorms_sh) * 0.05);
sensSet_sh     = find(sensUsed_sh);

bothSet_sh     = intersect(actSet_sh(:), sensSet_sh(:));
actOnlySet_sh  = setdiff(actSet_sh(:),  bothSet_sh);
sensOnlySet_sh = setdiff(sensSet_sh(:), bothSet_sh);
inactiveSet_sh = setdiff((1:nN_sh)', union(actSet_sh(:), sensSet_sh(:)));

scatter(nodePos_sh(inactiveSet_sh ,1), nodePos_sh(inactiveSet_sh ,2), 28, [0.75 0.75 0.75], 'filled', 'MarkerEdgeColor', 'none');
scatter(nodePos_sh(sensOnlySet_sh ,1), nodePos_sh(sensOnlySet_sh ,2), 42, cb_trunc,  'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.7);
scatter(nodePos_sh(actOnlySet_sh  ,1), nodePos_sh(actOnlySet_sh  ,2), 42, cb_diag,   'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.7);
scatter(nodePos_sh(bothSet_sh     ,1), nodePos_sh(bothSet_sh     ,2), 56, cb_purple, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.9);

title('(f) IEEE 13-bus EA controller', 'FontSize', titFS, 'Interpreter', 'latex');

exportgraphics(fig1, 'perf_bounds_convergence.pdf', 'ContentType', 'vector');
saveas(fig1, 'perf_bounds_convergence.png');
fprintf('Figure 1 saved → perf_bounds_convergence.*\n');

%% ===================== Figure 2: Unstable System (Baseline vs Gershgorin) =====================
gens_u = (1:U.maxGen)';
cb_u = [0.00 0.45 0.74;    % blue   — Baseline
        0.85 0.33 0.10];   % orange — Gershgorin
ca_u = [0.75 0.85 1.00;
        1.00 0.82 0.72];

fig2 = figure('Position', [60 60 1400 420], 'Color', 'w');
set(fig2, 'DefaultAxesFontSize', axFS, 'DefaultAxesFontName', 'Times New Roman');

% --- (a) Best Cost (normalised by Dense LQR total cost) ---
subplot(1, 3, 1); hold on;
for c = 1:U.nCond
    mu = mean(U.bestCostAll(:,:,c), 2) / U.J_dense_ref;
    sd = std(U.bestCostAll(:,:,c),  0, 2) / U.J_dense_ref;
    fill([gens_u; flipud(gens_u)], [mu+sd; flipud(mu-sd)], ...
        ca_u(c,:), 'EdgeColor', 'none', 'FaceAlpha', 0.45);
end
h = gobjects(U.nCond,1);
for c = 1:U.nCond
    h(c) = plot(gens_u, mean(U.bestCostAll(:,:,c),2) / U.J_dense_ref, ...
        '-', 'Color', cb_u(c,:), 'LineWidth', lw);
end
hRef = yline(1.0, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.8);
xlabel('Generation', 'FontSize', labFS);
ylabel('Best Cost $/ \; J_{\mathrm{dense}}$', 'FontSize', labFS, 'Interpreter', 'latex');
title('(a) Best Cost Convergence', 'FontSize', titFS);
legend([h; hRef], [U.condNames, {'Dense LQR ($=1$)'}], ...
    'Location', 'northeast', 'FontSize', legFS, 'Interpreter', 'latex');
grid on; box on;

% --- (b) Unstable Count ---
subplot(1, 3, 2); hold on;
for c = 1:U.nCond
    mu = mean(U.unstableAll(:,:,c), 2);
    sd = std(U.unstableAll(:,:,c),  0, 2);
    fill([gens_u; flipud(gens_u)], [mu+sd; flipud(mu-sd)], ...
        ca_u(c,:), 'EdgeColor', 'none', 'FaceAlpha', 0.45);
end
for c = 1:U.nCond
    plot(gens_u, mean(U.unstableAll(:,:,c),2), '-', 'Color', cb_u(c,:), 'LineWidth', lw);
end
xlabel('Generation', 'FontSize', labFS);
ylabel('Unstable Count', 'FontSize', labFS);
title('(b) Stability Failures', 'FontSize', titFS);
legend(U.condNames, 'Location', 'northeast', 'FontSize', legFS);
ylim([0 U.ea_params.popSize]);
grid on; box on;

% --- (c) Gershgorin Repairs ---
subplot(1, 3, 3); hold on;
mu_r = mean(U.repairAll(:,:,2), 2);
sd_r = std(U.repairAll(:,:,2),  0, 2);
fill([gens_u; flipud(gens_u)], [mu_r+sd_r; flipud(max(mu_r-sd_r,0))], ...
    ca_u(2,:), 'EdgeColor', 'none', 'FaceAlpha', 0.45);
plot(gens_u, mu_r, '-', 'Color', cb_u(2,:), 'LineWidth', lw);
xlabel('Generation', 'FontSize', labFS);
ylabel('Repairs / Gen', 'FontSize', labFS);
title('(c) Gershgorin Repairs', 'FontSize', titFS);
grid on; box on;

sgtitle(sprintf('EA-LQR: Baseline vs Gershgorin  (Grid %d\\times%d, \\rho(A)=%.2f, %d seeds)', ...
    U.gridSize, U.gridSize, U.rho_A, U.numSeeds), ...
    'FontSize', 16, 'FontWeight', 'bold', 'FontName', 'Times New Roman');

exportgraphics(fig2, 'perf_bounds_unstable.pdf', 'ContentType', 'vector');
saveas(fig2, 'perf_bounds_unstable.png');
fprintf('Figure 2 (unstable) saved → perf_bounds_unstable.*\n');
