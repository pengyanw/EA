clear
close all
addpath(genpath(pwd));

%%
gs = {@() simple_2d_graph(10), @() dc_network()};
names = {'2d-mesh', 'dc-118ieee'};
eta_press = {[-1:-1:-4], [-1:-1:-4]};
dts = [1, 0.005];
kappas = 0:1:5;

% Storage for combined Fig. 3
all_regrets = cell(1, 2);
all_errors  = cell(1, 2);

%%
for idx = 1:2
    g = gs{idx}();
    name = names{idx};
    dt = dts(idx);
    eta_pre = eta_press{idx};
    fprintf('%s\n', name);
    [n, ~] = size(g.Nodes);

    problem.name = names{idx};
    problem.graph = g;

    I = eye(n);
    if strcmp(name, '2d-mesh')
        L = laplacian_matrix(g);
        A = [I, dt*I; zeros(n), I - dt*L];
    elseif strcmp(name, 'dc-118ieee')
        L = laplacian_matrix(g);
        A = [I, dt*I; -dt*L, I];
    end
    R = eye(n);

    Is = arrayfun(@(i) [i, n+i], 1:n, 'UniformOutput', false);
    Js = num2cell((1:n)');

    eta_vals = 2.^eta_pre;
    regrets  = zeros(length(eta_vals), length(kappas));
    sradiuss = zeros(size(regrets));
    errors   = zeros(size(regrets));

    dt_str = strrep(num2str(dt), '.', 'p');

    % Prepare save directories
    save_dir_K = fullfile('trunc-K');
    save_dir_P = fullfile('trunc-P', name);
    if ~exist(save_dir_K, 'dir'); mkdir(save_dir_K); end
    if ~exist(save_dir_P, 'dir'); mkdir(save_dir_P); end

    for i = 1:length(eta_vals)
        eta      = eta_vals(i);
        eta_log2 = eta_pre(i);
        fprintf('eta = %.4f\n', eta);

        B = [zeros(n); dt * eta * eye(n)];
        C = [eta*eye(n), zeros(n)];
        Q = C' * C;

        % Solve discrete-time ARE
        [P, ~, ~] = idare(A, B, Q, R);
        K = (R + B' * P * B) \ (B' * P * A);

        % Fix 1: Save system matrices per eta (matches Julia matwrite)
        filename_sys = sprintf('%s-eta%d.mat', name, eta_log2);
        save(filename_sys, 'A', 'B', 'C', 'Q', 'P', 'K');

        for j = 1:length(kappas)
            kappa  = kappas(j);
            K_trunc = truncate(K, g, Is, Js, kappa, problem);

            Acl = A - B * K_trunc;
            s   = max(abs(eig(Acl)));
            sradiuss(i, j) = s;

            if s < 1
                P_trunc = dlyap(Acl', Q + K_trunc' * R * K_trunc);
                regrets(i, j) = norm(P_trunc - P, 2) / norm(P, 2);

                % Fix 3: save P_trunc (variable name matches Julia)
                filename_P = sprintf('P_eta%d_dt%s_kappa%d.mat', eta_log2, dt_str, kappa);
                save(fullfile(save_dir_P, filename_P), 'P_trunc');
            else
                regrets(i, j) = Inf;
            end

            errors(i, j) = norm(K - K_trunc, 2) / norm(K, 2);

            % Fix 2: save K_trunc to trunc-K/ (path and variable match Julia)
            filename_K = sprintf('K_eta%d_dt%s_kappa%d.mat', eta_log2, dt_str, kappa);
            save(fullfile(save_dir_K, filename_K), 'K_trunc');
        end
    end

    % Fix 4: Save aggregate regret_data and error_data (4×6, matches Julia)
    if ~exist('regret_data', 'dir'); mkdir('regret_data'); end
    if ~exist('error_data',  'dir'); mkdir('error_data');  end

    regret = regrets;  % 4×6
    save(fullfile('regret_data', [name, '.mat']), 'regret');

    % Save error with key name 'error' (matches Julia Dict key) without
    % shadowing MATLAB built-in by using struct save form
    s_err.error = errors;  % 4×6
    save(fullfile('error_data', [name, '.mat']), '-struct', 's_err');

    % Store for combined Fig. 3
    all_regrets{idx} = regrets;
    all_errors{idx}  = errors;

    % --- Individual PDFs (Fix 5a): matching Julia savefig naming ---
    markers = {'o', '+', 's', 'd'};

    % Individual P-gap figure
    fig_p = figure('Position', [100, 100, 450, 300]);
    for k = 1:length(eta_vals)
        semilogy(kappas, regrets(k,:), [markers{k}, '-'], ...
            'LineWidth', 1.5, 'MarkerSize', 7);
        hold on;
    end
    hold off;
    box on;
    xlabel('$\kappa$', 'Interpreter', 'latex', 'FontSize', 12);
    ylabel('$\|\mathbf{P}^\kappa - \mathbf{P}^\star\| / \|\mathbf{P}^\star\|$', ...
        'Interpreter', 'latex', 'FontSize', 12);
    legend(arrayfun(@(e) sprintf('$\\eta = 2^{%d}$', e), eta_pre, 'UniformOutput', false), ...
        'Interpreter', 'latex', 'Location', 'southwest', 'FontSize', 11);
    set(gca, 'FontSize', 12);
    exportgraphics(fig_p, [name, '-P.pdf'], 'ContentType', 'vector');

    % Individual K-error figure (errors' is 6×4)
    fig_k = figure('Position', [100, 100, 450, 300]);
    errors_plot = errors';  % 6×4
    for k = 1:length(eta_vals)
        semilogy(kappas, errors_plot(:,k), [markers{k}, '-'], ...
            'LineWidth', 1.5, 'MarkerSize', 7);
        hold on;
    end
    hold off;
    box on;
    xlabel('$\kappa$', 'Interpreter', 'latex', 'FontSize', 12);
    ylabel('$\|\mathbf{K}^\kappa - \mathbf{K}^\star\| / \|\mathbf{K}^\star\|$', ...
        'Interpreter', 'latex', 'FontSize', 12);
    legend(arrayfun(@(e) sprintf('$\\eta = 2^{%d}$', e), eta_pre, 'UniformOutput', false), ...
        'Interpreter', 'latex', 'Location', 'southwest', 'FontSize', 11);
    set(gca, 'FontSize', 12);
    exportgraphics(fig_k, [name, '-K.pdf'], 'ContentType', 'vector');
end

%% Fix 5b: Combined 2×2 Fig. 3 (paper figure)
markers       = {'o', '+', 's', 'd'};
panel_titles  = {'2d-mesh', 'dc-118ieee'};
eta_pre_all   = [-1, -2, -3, -4];
lbl_strings   = arrayfun(@(e) sprintf('$\\eta = 2^{%d}$', e), eta_pre_all, 'UniformOutput', false);

fig3 = figure('Position', [100, 100, 900, 600]);
for idx = 1:2
    reg  = all_regrets{idx};   % 4×6
    errs = all_errors{idx}';   % 6×4

    % Top row: P-gap
    subplot(2, 2, idx);
    for k = 1:4
        semilogy(kappas, reg(k,:), [markers{k}, '-'], ...
            'LineWidth', 1.5, 'MarkerSize', 7);
        hold on;
    end
    hold off;
    box on;
    title(panel_titles{idx}, 'FontSize', 12, 'Interpreter', 'none');
    xlabel('$\kappa$', 'Interpreter', 'latex', 'FontSize', 12);
    ylabel('$\|\mathbf{P}^\kappa\!-\!\mathbf{P}^\star\|\,/\,\|\mathbf{P}^\star\|$', ...
        'Interpreter', 'latex', 'FontSize', 11);
    legend(lbl_strings, 'Interpreter', 'latex', 'Location', 'southwest', 'FontSize', 9);
    set(gca, 'FontSize', 11);

    % Bottom row: K-error
    subplot(2, 2, idx + 2);
    for k = 1:4
        semilogy(kappas, errs(:,k), [markers{k}, '-'], ...
            'LineWidth', 1.5, 'MarkerSize', 7);
        hold on;
    end
    hold off;
    box on;
    title(panel_titles{idx}, 'FontSize', 12, 'Interpreter', 'none');
    xlabel('$\kappa$', 'Interpreter', 'latex', 'FontSize', 12);
    ylabel('$\|\mathbf{K}^\kappa\!-\!\mathbf{K}^\star\|\,/\,\|\mathbf{K}^\star\|$', ...
        'Interpreter', 'latex', 'FontSize', 11);
    legend(lbl_strings, 'Interpreter', 'latex', 'Location', 'southwest', 'FontSize', 9);
    set(gca, 'FontSize', 11);
end
exportgraphics(fig3, 'fig3.pdf', 'ContentType', 'vector');
