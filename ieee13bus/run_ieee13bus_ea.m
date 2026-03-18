%% run_ieee13bus_ea.m
% Apply EA-LQR with Gershgorin repair to the IEEE 13-bus test feeder.
% Compare Dense LQR, Diagonal (local) LQR, and EA-LQR.
%
% Two scenarios:
%   A) Full actuation  (all 13 buses)
%   B) Partial actuation (4 key DER buses)
%
% Author: Pengyang Wu
% Date: 2026-02-17

clear; clc; close all;
addpath(genpath(pwd));

%% ===================== Configuration =====================
eaSeeds  = [1, 17, 42];
numSeeds = numel(eaSeeds);

ea_params.popSize  = 20;
ea_params.maxGen   = 150;
ea_params.pMutate  = 0.05;
ea_params.pCross   = 0.8;
ea_params.nTop     = 10;
ea_params.alpha    = 0;
maxGen = ea_params.maxGen;

opts.verbose       = false;
opts.useGersRepair = true;

% Two actuation scenarios
scenarios = {
    struct('name', 'Full (13 buses)',   'buses', 1:13);
    %struct('name', 'Partial (4 DERs)',  'buses', [1, 4, 7, 13]);  % 650, 634, 671, 675
};
nScen = numel(scenarios);

axFS = 13;  labFS = 14;  titFS = 15;  legFS = 11;  lw = 2.5;
cb_ea   = [0.00 0.45 0.74];
cb_dens = [0.50 0.50 0.50];
cb_diag = [0.85 0.33 0.10];
ca_ea   = [0.75 0.85 1.00];

%% ===================== Figure 0: Topology =====================
cfg0 = struct('actuatedBuses', 1:13);
[~, binfo] = build_ieee13bus_system(cfg0);

fig0 = figure('Position', [60 500 480 400], 'Color', 'w');
hold on;
coords = binfo.nodeCoords;
ed = binfo.edgeData;
for e = 1:size(ed, 1)
    i = ed(e,1);  j = ed(e,2);
    plot([coords(i,1) coords(j,1)], [coords(i,2) coords(j,2)], ...
        'k-', 'LineWidth', 1.5);
end
scatter(coords(:,1), coords(:,2), 120, cb_ea, 'filled', 'MarkerEdgeColor', 'k');
for b = 1:binfo.nBus
    text(coords(b,1)+0.15, coords(b,2)+0.2, binfo.busNames{b}, ...
        'FontSize', 10, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
end


%% ===================== Main Loop =====================
% Storage
convCurves = cell(nScen, numSeeds);
J_dense    = zeros(nScen, 1);
J_diag     = zeros(nScen, 1);
J_ea       = zeros(nScen, numSeeds);
nnz_dense  = zeros(nScen, 1);
nnz_diag   = zeros(nScen, 1);
nnz_ea     = zeros(nScen, numSeeds);
rho_A_vec  = zeros(nScen, 1);
nxVec      = zeros(nScen, 1);
nuVec      = zeros(nScen, 1);

for sc = 1:nScen
    fprintf('\n========== Scenario %d: %s ==========\n', sc, scenarios{sc}.name);

    cfg.Ts = 0.2;
    cfg.actuatedBuses = scenarios{sc}.buses;
    [sys, binfo] = build_ieee13bus_system(cfg);
    A  = sys.A;   B = sys.B2;
    Nx = sys.Nx;  Nu = sys.Nu;
    nxVec(sc) = Nx;  nuVec(sc) = Nu;
    max_links = Nu * Nx;
    rho_A_vec(sc) = max(abs(eig(A)));

    fprintf('  Nx=%d, Nu=%d, rho(A)=%.4f\n', Nx, Nu, rho_A_vec(sc));

    Q_bm = eye(Nx);  R_bm = eye(Nu);
    K_full = -dlqr(A, B, Q_bm, R_bm);
    costBM = get_lqr_cost(A, B, Q_bm, R_bm, K_full);
    fprintf('  Benchmark (full dlqr) cost: %.4f\n', costBM);

    % --- Dense LQR (|K| >= 1e-3) ---
    K_dens = K_full;
    %K_dens(abs(K_dens) < 1e-3) = 0;
    J_dense(sc) = cost_EA(A, B, Q_bm, R_bm, K_dens, costBM, ea_params.alpha, max_links);
    nnz_dense(sc) = nnz(K_dens);
    fprintf('  Dense LQR: J=%.2f, nnz=%d/%d (%.1f%%)\n', ...
        J_dense(sc), nnz_dense(sc), max_links, 100*nnz_dense(sc)/max_links);

    % --- Diagonal LQR: actuator k -> states [2*bus-1, 2*bus] ---
    K_dig = zeros(Nu, Nx);
    actBuses = scenarios{sc}.buses;
    for k = 1:Nu
        bIdx = actBuses(k);
        c1 = 2*bIdx - 1;  c2 = 2*bIdx;
        if c2 <= Nx
            K_dig(k, c1) = K_full(k, c1);
            K_dig(k, c2) = K_full(k, c2);
        elseif c1 <= Nx
            K_dig(k, c1) = K_full(k, c1);
        end
    end
    if max(abs(eig(A + B * K_dig))) >= 1.0
        [K_dig, info_d] = gersgorin_stabilize_K(A, B, K_dig, K_full);
        if ~info_d.success
            fprintf('  WARNING: Diagonal K still unstable\n');
        end
    end
    J_diag(sc) = cost_EA(A, B, Q_bm, R_bm, K_dig, costBM, ea_params.alpha, max_links);
    nnz_diag(sc) = nnz(K_dig);
    fprintf('  Diagonal LQR: J=%.2f, nnz=%d\n', J_diag(sc), nnz_diag(sc));

    % --- EA runs ---
    for s = 1:numSeeds
        seed = eaSeeds(s);
        rng(seed);
        fprintf('  EA seed %d ...', seed);

        [result, history] = ea_lqr_codesign_gershgorin(A, B, ea_params, opts);
        convCurves{sc, s} = history.bestCost;
        J_ea(sc, s)   = result.cost;
        nnz_ea(sc, s) = nnz(result.K_sparse);

        fprintf(' J=%.2f, nnz=%d, act=%d/%d, sens=%d/%d\n', ...
            J_ea(sc,s), nnz_ea(sc,s), ...
            numel(result.active_actuators), Nu, ...
            numel(result.active_sensors), Nx);
    end
end

%% ===================== Console Summary =====================
fprintf('\n================ IEEE 13-Bus Summary ================\n');
for sc = 1:nScen
    Jref = J_dense(sc);
    fprintf('\nScenario: %s (Nx=%d, Nu=%d, rho(A)=%.3f)\n', ...
        scenarios{sc}.name, nxVec(sc), nuVec(sc), rho_A_vec(sc));
    fprintf('  %-18s  J/J_dense  nnz\n', 'Method');
    fprintf('  %s\n', repmat('-', 1, 42));
    fprintf('  %-18s  %7.3f    %d\n', 'Dense LQR', 1.000, nnz_dense(sc));
    if J_diag(sc) < 1e8
        fprintf('  %-18s  %7.3f    %d\n', 'Diagonal LQR', J_diag(sc)/Jref, nnz_diag(sc));
    else
        fprintf('  %-18s  %7s    %d\n', 'Diagonal LQR', 'UNSTAB', nnz_diag(sc));
    end
    for s = 1:numSeeds
        fprintf('  %-18s  %7.3f    %d\n', ...
            sprintf('EA (seed %d)', eaSeeds(s)), J_ea(sc,s)/Jref, nnz_ea(sc,s));
    end
    fprintf('  %-18s  %7.3f +/- %.3f\n', 'EA (mean)', ...
        mean(J_ea(sc,:))/Jref, std(J_ea(sc,:))/Jref);

    imp_dense = (J_dense(sc) - mean(J_ea(sc,:))) / J_dense(sc) * 100;
    if J_diag(sc) < 1e8
        imp_diag = (J_diag(sc) - mean(J_ea(sc,:))) / J_diag(sc) * 100;
    else
        imp_diag = NaN;
    end
    fprintf('  EA improvement: vs Dense %+.1f%%, vs Diagonal %+.1f%%\n', imp_dense, imp_diag);
end

%% ===================== Figure 1: Convergence =====================
fig1 = figure('Position', [60 60 1100 420], 'Color', 'w');
set(fig1, 'DefaultAxesFontSize', axFS, 'DefaultAxesFontName', 'Times New Roman');
gens = (1:maxGen)';

for sc = 1:nScen
    subplot(1, nScen, sc); hold on;
    Jref = J_dense(sc);

    % EA curves normalized by Dense
    normMat = zeros(maxGen, numSeeds);
    for s = 1:numSeeds
        normMat(:, s) = convCurves{sc, s} / Jref;
    end
    mu = mean(normMat, 2);
    sd = std(normMat, 0, 2);

    fill([gens; flipud(gens)], [mu+sd; flipud(max(mu-sd, 0))], ...
        ca_ea, 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    hEA = plot(gens, mu, '-', 'Color', cb_ea, 'LineWidth', lw);

    hDens = yline(1.0, '--', 'Color', cb_dens, 'LineWidth', 2);

    J_diag_n = J_diag(sc) / Jref;
    if J_diag(sc) < 1e8
        hDiag = yline(J_diag_n, '-.', 'Color', cb_diag, 'LineWidth', 2);
    else
        hDiag = plot(NaN, NaN, '-.', 'Color', cb_diag, 'LineWidth', 2);
        text(maxGen/2, 1.4, 'Diagonal: Unstable', 'Color', cb_diag, ...
            'FontSize', legFS, 'FontName', 'Times New Roman', ...
            'HorizontalAlignment', 'center');
    end

    % Mark gen where EA < Dense
    idx_cross = find(mu < 1.0, 1, 'first');
    if ~isempty(idx_cross)
        plot(idx_cross, mu(idx_cross), 'v', 'MarkerSize', 8, ...
            'MarkerFaceColor', cb_ea, 'MarkerEdgeColor', 'k');
    end

    xlabel('Generation', 'FontSize', labFS);
    if sc == 1
        ylabel('$J \;/\; J_{\mathrm{dense}}$', 'FontSize', labFS, 'Interpreter', 'latex');
    end
    title(sprintf('(%c) %s', char('a'+sc-1), scenarios{sc}.name), ...
        'FontSize', titFS, 'FontName', 'Times New Roman');
    grid on; box on;

    if sc == nScen
        legend([hEA, hDens, hDiag], ...
            {'EA-LQR', 'Dense LQR ($=1$)', 'Diagonal LQR'}, ...
            'Location', 'northeast', 'FontSize', legFS, 'Interpreter', 'latex');
    end
end

sgtitle('IEEE 13-Bus: EA-LQR Convergence (/ Dense LQR)', ...
    'FontSize', 16, 'FontWeight', 'bold', 'FontName', 'Times New Roman');

%% ===================== Figure 2: Sparsity Pattern =====================
% Show K_sparse structure for the last EA run of each scenario
fig2 = figure('Position', [60 520 1100 380], 'Color', 'w');
set(fig2, 'DefaultAxesFontSize', axFS, 'DefaultAxesFontName', 'Times New Roman');

for sc = 1:nScen
    % Re-run with last seed to get K
    cfg.Ts = 0.2;
    cfg.actuatedBuses = scenarios{sc}.buses;
    [sys, ~] = build_ieee13bus_system(cfg);
    rng(eaSeeds(end));
    [result, ~] = ea_lqr_codesign_gershgorin(sys.A, sys.B2, ea_params, opts);

    subplot(1, nScen, sc);
    imagesc(abs(result.K_sparse));
    colormap(flipud(hot));
    colorbar('FontSize', axFS-2);
    xlabel('State index', 'FontSize', labFS);
    ylabel('Actuator index', 'FontSize', labFS);
    title(sprintf('(%c) %s — nnz=%d', char('a'+sc-1), scenarios{sc}.name, ...
        nnz(result.K_sparse)), 'FontSize', titFS, 'FontName', 'Times New Roman');
    axis tight;
end

sgtitle('EA-LQR Sparse Controller Structure on IEEE 13-Bus', ...
    'FontSize', 16, 'FontWeight', 'bold', 'FontName', 'Times New Roman');

%% ===================== Save =====================
outDir = fullfile(fileparts(mfilename('fullpath')), 'results');
if ~exist(outDir, 'dir'), mkdir(outDir); end

saveas(fig0, fullfile(outDir, 'ieee13bus_topology.png'));
exportgraphics(fig1, fullfile(outDir, 'ieee13bus_convergence.pdf'), 'ContentType', 'vector');
saveas(fig1, fullfile(outDir, 'ieee13bus_convergence.png'));
exportgraphics(fig2, fullfile(outDir, 'ieee13bus_sparsity.pdf'), 'ContentType', 'vector');
saveas(fig2, fullfile(outDir, 'ieee13bus_sparsity.png'));

save(fullfile(outDir, 'ieee13bus_results.mat'), ...
    'scenarios', 'eaSeeds', 'ea_params', 'convCurves', ...
    'J_dense', 'J_diag', 'J_ea', 'nnz_dense', 'nnz_diag', 'nnz_ea', 'rho_A_vec');

fprintf('\nResults saved to %s\n', outDir);
fprintf('========== IEEE 13-Bus Experiment Complete ==========\n');

