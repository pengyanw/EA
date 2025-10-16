%% Random Forest Regression from EA popbuffer
clear; clc;
addpath(genpath(pwd));  
close all

% --------------------------------------------------
% 1. 加载数据
% --------------------------------------------------
load("G:\其他计算机\God desk\Research\Lisali\EA_new_and_old\EA\cache\popbuffer.mat");   % 若 popbuffer 已在 workspace，可省略
[nPop, nGen] = size(popbuffer);
% fprintf('Loaded popbuffer: %d populations × %d generations\n', nPop, nGen);
% %% 1. System Generation (No changes here)
% % =========================================================================
% gridSize       = 5;
% connectThresh  = 0.5;
% Ts             = 0.2;
% actDensity     = 1;
% seed           = 17; 
% 
% numNodes    = gridSize*gridSize;
% [adjMtx, nodeCoords, susceptMtx, inertiasInv, dampings] = generate_grid_topology(gridSize, connectThresh, seed);
% plot_graph(adjMtx, nodeCoords, 'k');
% 
% numActs       = round(actDensity*numNodes);
% actuatedNodes = randsample(numNodes, numActs);    
% sys           = generate_grid_plant(actuatedNodes, adjMtx, susceptMtx, inertiasInv, dampings, Ts);
% 
% Nx = sys.Nx;
% Nu = sys.Nu;
% A  = sys.A;
% B_ = sys.B2;


%% 

% --------------------------------------------------
% 2. 展开所有样本
% --------------------------------------------------
allData = cell2mat(reshape(popbuffer, nPop*nGen, []));   % [3000 × 76]

% --------------------------------------------------
% 3. 构造特征与输出
% --------------------------------------------------
X = allData(:, 1:75);    % 特征（Q+R）
Y = allData(:, 76);      % 输出（mask数量）

fprintf('Dataset: %d samples, %d features\n', size(X,1), size(X,2));

% --------------------------------------------------
% 4. 划分训练与测试集
% --------------------------------------------------
cv = cvpartition(size(X,1), 'HoldOut', 0.2);
X_train = X(training(cv), :);
Y_train = Y(training(cv), :);
X_test  = X(test(cv),:);
Y_test  = Y(test(cv), :);

% --------------------------------------------------
% 5. 随机森林训练
% --------------------------------------------------
numTrees = 100;   % 可调参数
rfModel = TreeBagger(numTrees, X_train, Y_train, ...
    'Method', 'regression', ...
    'MinLeafSize', 5, ...
    'OOBPrediction', 'on', ...
    'OOBPredictorImportance', 'on');
alpha = 0.5;
% for i = 1:numTrees
%     ind = randperm(75, 75*alpha)
%     features(i,:) = ind
%     Xtrain = Xtrain(:, ind)
%     tree = Decisiontree(Xtrain, ytrain, Minleafsize=5)
% end

% --------------------------------------------------
% 6. 评估模型性能
% --------------------------------------------------
Y_pred = predict(rfModel, X_test);
mse_val = mean((Y_test - Y_pred).^2);
fprintf('Test MSE = %.4f\n', mse_val);

% --------------------------------------------------
% 7. 特征重要性分析
% --------------------------------------------------
figure;
bar(rfModel.OOBPermutedPredictorDeltaError);
xlabel('Feature Index (1–75)');
ylabel('Importance');
title('Random Forest Feature Importance');
grid on;
if ~exist('figures1', 'dir')
    mkdir('figures1');
end


% --------------------------------------------------
% 8. 保存模型与结果
% --------------------------------------------------
if ~exist('cache', 'dir')
    mkdir('cache');
end
save('cache/rf_model.mat', 'rfModel', 'mse_val');

% --------------------------------------------------
% 9. 可选：查看 OOB 误差曲线
% --------------------------------------------------
figure;
oobErrorBaggedEnsemble = oobError(rfModel);
plot(oobErrorBaggedEnsemble, 'LineWidth', 1.5);
xlabel('Number of grown trees');
ylabel('Out-of-bag MSE');
title('OOB Error Progress');
grid on;
% 10. 保存图像到 figures1 文件夹
% --------------------------------------------------
saveas(figure(1), fullfile('figures1', 'rf_feature_importance.png'));
saveas(figure(2), fullfile('figures1', 'rf_oob_error.png'));

fprintf('All figures saved to ./figures1/\n');