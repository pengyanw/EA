%% Random Forest training for K importance (simplified features)
clear; clc; close all;
addpath(genpath(pwd));

%% ========================= 1. Parameters ================================
numSamples   = 1000;       % number of training systems
gridSize     = 5;         % grid size
connectThresh= 0.5;
Ts           = 0.2;
actDensity   = 1.0;
QscaleRange  = [0.1, 10];
RscaleRange  = [0.1, 10];

fprintf('Generating %d training samples...\n', numSamples);

%% ========================= 2. Data Buffers ==============================
X_all = [];
Y_all = [];

%% ========================= 3. Generate Random Systems ===================
for k = 1:numSamples
    rng(k);
    % --- generate topology & plant
    [adjMtx,~, susceptMtx, inertiasInv, dampings] = ...
        generate_grid_topology(gridSize, connectThresh, k);
    numNodes = gridSize^2;
    numActs  = round(actDensity * numNodes);
    actNodes = randsample(numNodes, numActs);
    sys = generate_grid_plant(actNodes, adjMtx, susceptMtx, inertiasInv, dampings, Ts);
    A = sys.A; B = sys.B2;
    Nx = sys.Nx; Nu = sys.Nu;

    % --- random Q,R
    diagQ = QscaleRange(1) + diff(QscaleRange)*rand(1,Nx);
    diagR = RscaleRange(1) + diff(RscaleRange)*rand(1,Nu);
    Q = diag(diagQ);
    R = diag(diagR);

    % --- dense LQR controller
    try
        K_dense = -dlqr(A,B,Q,R);
    catch
        continue;
    end

        % ---------- 预计算图距离 ----------
    G = graph(adjMtx);
    D = distances(G);
    D(isinf(D)) = max(D(~isinf(D))) + 1;
    D_expanded = kron(D, ones(1, Nx/numNodes));   % 维度匹配到 Nu x Nx
    
    % ---------- 全局尺度 ----------
    normA = norm(A, 'fro');
    normB = norm(B, 'fro');
    
    % ---------- AB 的局部特征 ----------
    colA  = sum(abs(A), 1)';       % Nx x 1   每个状态列的1范数
    diagA = abs(diag(A));          % Nx x 1   对角项
    bPow  = vecnorm(B, 2, 1)';     % Nu x 1   每个输入通道的2范数
    
    Galign = abs(B' * A);          % Nu x Nx  输入与第j列动力学的对齐度
    
    % ---------- 展开到 (i,j) 配对 ----------
    [II, JJ] = ndgrid(1:Nu, 1:Nx);  % II: actuator index, JJ: state index
    
    feat_d      = D_expanded(:);
    
    feat_colA   = colA(JJ(:));      % 对应 j
    feat_diagA  = diagA(JJ(:));     % 对应 j
    feat_bPow   = bPow(II(:));      % 对应 i
    feat_align  = Galign(:);        % 对应 (i,j)
    
    % ---------- 拼装最终特征矩阵 ----------
    X = [feat_d, feat_colA, feat_diagA, feat_bPow, feat_align];
    
    % 标签仍然用 |K_dense(i,j)|
    Y = abs(K_dense(:));


    X_all = [X_all; X];
    Y_all = [Y_all; Y];
end

fprintf('Total pairs collected: %d\n', size(X_all,1));

%% ========================= 4. Train Random Forest =======================
numTrees = 100;
fprintf('Training Random Forest Regressor...\n');
rfModel = TreeBagger(numTrees, X_all, Y_all, ...
    'Method','regression', ...
    'OOBPrediction','on', ...
    'OOBPredictorImportance','on');

%% ========================= 5. Evaluate Model ============================
Y_pred = oobPredict(rfModel);
mse_val = mean((Y_all - Y_pred).^2);
fprintf('OOB MSE = %.4e\n', mse_val);

figure;
bar(rfModel.OOBPermutedPredictorDeltaError);
xticklabels({'d(i,j)','||A||_F','||B||_F'});
xlabel('Feature');
ylabel('Importance');
title('Random Forest Feature Importance (Simplified)');
grid on;

%% ========================= 6. Save Model ================================
if ~exist('cache','dir'), mkdir('cache'); end
save('cache/rf_Kimportance_simple.mat','rfModel');

%% ========================= 7. Predict on new graph ======================
fprintf('Testing on a new graph...\n');
seed = 5; rng(seed);
[adjMtx, nodeCoords, susceptMtx, inertiasInv, dampings] = ...
    generate_grid_topology(gridSize, connectThresh, seed);
numNodes = gridSize^2;
actNodes = randsample(numNodes, numActs);
sys = generate_grid_plant(actNodes, adjMtx, susceptMtx, inertiasInv, dampings, Ts);
A = sys.A; B = sys.B2; Nx = sys.Nx; Nu = sys.Nu;

G = graph(adjMtx);
D = distances(G);
D(isinf(D)) = max(D(~isinf(D)))+1e2;
D_expanded = kron(D, ones(1, Nx/numNodes));
feat_d      = D_expanded(:);
feat_colA   = colA(JJ(:));      % 对应 j
feat_diagA  = diagA(JJ(:));     % 对应 j
feat_bPow   = bPow(II(:));      % 对应 i
feat_align  = Galign(:);        % 对应 (i,j)

% ---------- 拼装最终特征矩阵 ----------
X_new = [feat_d, feat_colA, feat_diagA, feat_bPow, feat_align];
    
Y_pred_new = predict(rfModel, X_new);
K_pred = reshape(Y_pred_new, [Nu, Nx]);

figure;
imagesc(K_pred); colorbar;
xlabel('State index'); ylabel('Actuator index');
title('Predicted |K| Importance by RF (Simplified features)');
legend('feat_d', 'feat_colA', 'feat_diagA', 'feat_bPow', 'feat_align')
%% ============================================================
% 9. Validation across random seeds (RF vs Benchmark)
% =============================================================
nSeeds = 10;                 % number of random test systems
cost_RF  = zeros(nSeeds,1);
cost_BM  = zeros(nSeeds,1);

fprintf('\n==== Random-Forest Validation ====\n');

for s = 1:nSeeds
    % ---------- Generate new random system ----------
    gridSize = 5;
    connectThresh = 0.5;
    Ts = 0.2;
    actDensity = 1;
    seed = s * 2;           % different seed
    numNodes = gridSize^2;

    [adjMtx, nodeCoords, susceptMtx, inertiasInv, dampings] = ...
        generate_grid_topology(gridSize, connectThresh, seed);
    numActs       = round(actDensity*numNodes);
    actuatedNodes = randsample(numNodes, numActs);
    sys = generate_grid_plant(actuatedNodes, adjMtx, susceptMtx, inertiasInv, dampings, Ts);

    A = sys.A;
    B = sys.B2;
    Nx = sys.Nx; Nu = sys.Nu;

    % ---------- Benchmark dense LQR ----------
    Q_bm = eye(Nx);
    R_bm = eye(Nu);
    K_dense = -dlqr(A,B,Q_bm,R_bm);
    K_dense(abs(K_dense)<1e-3) = 0;
    costBM = get_lqr_cost(A, B, Q_bm, R_bm, K_dense);
    cost_BM(s) = cost_EA(A,B,Q_bm,R_bm,K_dense,costBM,0,Nu*Nx);

    % ---------- Construct features for RF ----------
    G = graph(adjMtx);
    D = distances(G);
    D(isinf(D)) = max(D(~isinf(D))) + 1;
  
    colA  = sum(abs(A),1)';
    diagA = abs(diag(A));
    bPow  = vecnorm(B,2,1)';
    Galign = abs(B'*A);
    [II,JJ] = ndgrid(1:Nu,1:Nx);

    feat_d      = D_expanded(:);

    feat_colA   = colA(JJ(:));      % 对应 j
    feat_diagA  = diagA(JJ(:));     % 对应 j
    feat_bPow   = bPow(II(:));      % 对应 i
    feat_align  = Galign(:);        % 对应 (i,j)

    X_test = [feat_d, feat_colA, feat_diagA, feat_bPow, feat_align];

    % ---------- Predict K and evaluate cost ----------
    Y_pred = predict(rfModel, X_test);
    K_pred = reshape(Y_pred, [Nu, Nx]);
    K_pred(abs(K_pred)<1e-3) = 0;
    cost_RF(s) = cost_EA(A,B,Q_bm,R_bm,K_pred,1,0,Nu*Nx);

    fprintf('Seed %2d:  RF cost = %.3f   |  Benchmark cost = %.3f\n', ...
             s, cost_RF(s), cost_BM(s));
end

% ----------  Plot comparison ----------
figure;
bar(1:nSeeds, [cost_BM, cost_RF]);
legend({'Benchmark LQR','Predicted K (RF)'},'Location','best');
xlabel('System Seed Index');
ylabel('Normalized Cost');
title('Cost Comparison: Random-Forest Prediction vs Benchmark');
grid on;

figure;
plot(1:nSeeds, cost_RF./cost_BM,'-o','LineWidth',1.5);
xlabel('System Seed Index');
ylabel('Cost Ratio (RF / Benchmark)');
title('Relative Performance of Random-Forest Predicted Gains');
grid on;

fprintf('\nAverage cost ratio (RF/Benchmark) = %.3f\n', mean(cost_RF./cost_BM));

%% ========================= 8. Visualization with log color scale =======
% Assume K_pred (from RF) and K_dense (from LQR) are already in workspace

% --- Normalize and clip for numerical safety
K_pred_abs  = abs(K_pred);
K_dense_abs = abs(K_dense);
K_pred_abs(K_pred_abs<=1e-3)   = 1e-6;
K_dense_abs(K_dense_abs==0) = 1e-6;

% --- Compute log-scale maps
K_pred_log  = log10(K_pred_abs);
K_dense_log = log10(K_dense_abs);

figure;
imagesc(K_pred_log);
colorbar;
xlabel('State index'); ylabel('Actuator index');
title('Predicted |K| (log_{10} scale)');
colormap(jet);
caxis([min(K_pred_log(:)) max(K_pred_log(:))]);
set(gca,'YDir','normal');
grid on;

% --- Compute error map
K_error = K_pred - K_dense;
K_error_abs = abs(K_error);
K_error_rel = K_error_abs ./ (abs(K_dense)+1e-6);

figure;
subplot(1,2,1);
imagesc(K_error);
colorbar;
xlabel('State index'); ylabel('Actuator index');
title('Signed Error (K_{pred} - K_{dense})');
colormap(jet); % if you have it, or jet otherwise
set(gca,'YDir','normal');

subplot(1,2,2);
imagesc(K_error_rel);
colorbar;
xlabel('State index'); ylabel('Actuator index');
title('Relative Error |ΔK| / |K_{dense}|');
colormap(jet);
set(gca,'YDir','normal');

% --- Summary metrics
mse_K = mean((K_pred(:)-K_dense(:)).^2);
mae_K = mean(abs(K_pred(:)-K_dense(:)));
fprintf('MSE between predicted and dense K = %.4e\n', mse_K);
fprintf('MAE between predicted and dense K = %.4e\n', mae_K);

%% ========================= 8. Visualization with log color scale =======
% Assume K_pred (from RF) and K_dense (from LQR) are already in workspace

% --- Normalize and clip for numerical safety
K_pred_abs  = abs(K_pred);
K_dense_abs = abs(K_dense);
K_pred_abs(K_pred_abs<=1e-3)   = 1e-6;
K_dense_abs(K_dense_abs==0)    = 1e-6;

% --- Compute log-scale maps
K_pred_log  = log10(K_pred_abs);
K_dense_log = log10(K_dense_abs);

figure(1);
imagesc(K_pred_log);
colorbar;
xlabel('State index'); ylabel('Actuator index');
title('Predicted |K| (log_{10} scale)');
colormap(jet);
caxis([min(K_pred_log(:)) max(K_pred_log(:))]);
set(gca,'YDir','normal');
grid on;

% --- Compute error map
K_error = K_pred - K_dense;
K_error_abs = abs(K_error);
K_error_rel = K_error_abs ./ (abs(K_dense)+1e-6);

figure(2);
subplot(1,2,1);
imagesc(K_error);
colorbar;
xlabel('State index'); ylabel('Actuator index');
title('Signed Error (K_{pred} - K_{dense})');
colormap(jet);
set(gca,'YDir','normal');

subplot(1,2,2);
imagesc(K_error_rel);
colorbar;
xlabel('State index'); ylabel('Actuator index');
title('Relative Error |ΔK| / |K_{dense}|');
colormap(jet);
set(gca,'YDir','normal');

% --- Summary metrics
mse_K = mean((K_pred(:)-K_dense(:)).^2);
mae_K = mean(abs(K_pred(:)-K_dense(:)));
fprintf('MSE between predicted and dense K = %.4e\n', mse_K);
fprintf('MAE between predicted and dense K = %.4e\n', mae_K);

%% ========================= 9. Save figures to designated folder ========
saveDir = 'E:\Research\Lisali\EA_new_and_old\EA\Reports\week8 figures';

% ensure the folder exists
if ~exist(saveDir, 'dir')
    mkdir(saveDir);
end

% timestamp for unique filenames
timestamp = datestr(now, 'yyyymmdd_HHMMSS');

% save log-scale predicted K
saveas(figure(1), fullfile(saveDir, ['K_pred_log_', timestamp, '.png']));
saveas(figure(1), fullfile(saveDir, ['K_pred_log_', timestamp, '.fig']));

% save error maps
saveas(figure(2), fullfile(saveDir, ['K_error_maps_', timestamp, '.png']));
saveas(figure(2), fullfile(saveDir, ['K_error_maps_', timestamp, '.fig']));

fprintf('Figures saved to: %s\n', saveDir);
