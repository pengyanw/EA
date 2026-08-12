%% plot_fig2_pdf.m
% Rebuild Figure 2 (unstable plant: Baseline vs Gershgorin repair) from saved
% data and export it as a vector PDF for the paper. No EA run required.
%
% Prerequisite:
%   run_multi_seed_unstable.m  ->  plot_only/fig_data_unstable.mat
%
% Usage:
%   cd plot_only
%   plot_fig2_pdf
%
% Output:
%   plot_only/perf_bounds_unstable.pdf   (vector, text selectable)
%   plot_only/perf_bounds_unstable.png
%
% Author: Pengyang Wu

clear; clc; close all;

matFile = 'fig_data_unstable.mat';
if ~exist(matFile, 'file')
    error('plot_fig2_pdf:missingData', ...
        ['%s not found. Run run_multi_seed_unstable.m first ' ...
         '(it writes this file into plot_only/).'], matFile);
end
U = load(matFile);

%% ===================== Style =====================
axFS = 18;  labFS = 19;  titFS = 19;  legFS = 15;  lw = 3;
cb = [0.00 0.45 0.74;    % blue   - Baseline (no repair)
      0.85 0.33 0.10];   % orange - Gershgorin repair
ca = [0.75 0.85 1.00;
      1.00 0.82 0.72];

gens = (1:U.maxGen)';

% Two panels, sized for a two-column paper figure. tiledlayout with compact
% spacing and tight padding leaves far less dead space than subplot, and
% exporting the layout object (rather than the figure) crops to the axes.
fig2 = figure('Position', [60 60 920 380], 'Color', 'w');
set(fig2, 'DefaultAxesFontSize', axFS, 'DefaultAxesFontName', 'Times New Roman');
tl = tiledlayout(fig2, 1, 2, 'TileSpacing', 'compact', 'Padding', 'tight');

%% --- (a) Best cost, normalized by the dense-LQR total cost ---
nexttile(tl); hold on;
for c = 1:U.nCond
    mu = mean(U.bestCostAll(:,:,c), 2) / U.J_dense_ref;
    sd = std(U.bestCostAll(:,:,c), 0, 2) / U.J_dense_ref;
    fill([gens; flipud(gens)], [mu+sd; flipud(max(mu-sd, 0))], ...
        ca(c,:), 'EdgeColor', 'none', 'FaceAlpha', 0.45);
end
for c = 1:U.nCond
    plot(gens, mean(U.bestCostAll(:,:,c), 2) / U.J_dense_ref, ...
        '-', 'Color', cb(c,:), 'LineWidth', lw);
end
yline(1.0, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.8);
xlabel('Generation', 'FontSize', labFS);
ylabel('Best Cost $/\;J_{\mathrm{dense}}$', 'FontSize', labFS, 'Interpreter', 'latex');
title('(a) Best Cost Convergence', 'FontSize', titFS);
% No legend here: panel (b) carries the same two conditions in the same colors.
% The dashed grey line at 1 is the dense-LQR reference -- state that in the
% caption, since it no longer has a legend entry.
grid on; box on;

%% --- (b) Unstable individuals per generation ---
nexttile(tl); hold on;
for c = 1:U.nCond
    mu = mean(U.unstableAll(:,:,c), 2);
    sd = std(U.unstableAll(:,:,c), 0, 2);
    fill([gens; flipud(gens)], [mu+sd; flipud(max(mu-sd, 0))], ...
        ca(c,:), 'EdgeColor', 'none', 'FaceAlpha', 0.45);
end
hb = gobjects(U.nCond, 1);
for c = 1:U.nCond
    hb(c) = plot(gens, mean(U.unstableAll(:,:,c), 2), '-', ...
        'Color', cb(c,:), 'LineWidth', lw);
end
xlabel('Generation', 'FontSize', labFS);
ylabel('Unstable Count', 'FontSize', labFS);
title('(b) Stability Failures', 'FontSize', titFS);
% Attach the legend to the line handles: the two fill() patches are drawn first,
% so a bare legend(condNames) would label the shaded bands instead.
legend(hb, U.condNames, 'Location', 'northeast', 'FontSize', legFS);
ylim([0 U.ea_params.popSize]);
grid on; box on;

% Panel (c) (repairs per generation) is intentionally omitted -- the repair
% counts are reported in the text instead. U.repairAll is still in the .mat if
% it is ever wanted back.

title(tl, sprintf(['EA-LQR: Baseline vs Gershgorin  ' ...
    '(Grid %d\\times%d, \\rho(A)=%.2f, %d seeds)'], ...
    U.gridSize, U.gridSize, U.rho_A, U.numSeeds), ...
    'FontSize', 16, 'FontWeight', 'bold', 'FontName', 'Times New Roman');

%% ===================== Export =====================
% ContentType 'vector' keeps text selectable and lines resolution-independent.
% Exporting the tiledlayout rather than the figure crops to the drawn content,
% removing the figure-window margin entirely.
exportgraphics(tl, 'perf_bounds_unstable.pdf', ...
    'ContentType', 'vector', 'BackgroundColor', 'white');
exportgraphics(tl, 'perf_bounds_unstable.png', 'Resolution', 200);
fprintf('Figure 2 saved -> perf_bounds_unstable.pdf (vector) and .png\n');
fprintf('  source: %s  (%d seeds, grid %dx%d, rho(A)=%.2f)\n', ...
    matFile, U.numSeeds, U.gridSize, U.gridSize, U.rho_A);
