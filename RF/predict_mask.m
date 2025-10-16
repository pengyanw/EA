

% 输入: N×75 的新样本矩阵 QR_input
% 输出: N×1 的预测 mask 数值
%
% 调用方式示例：
%   load('cache/rf_model.mat');
%   newQR = rand(5, 75);          % 示例新样本
%   y_pred = predict_mask(rfModel, newQR);
%
% 函数定义如下：
function y_pred = predict_mask(rfModel, QR_input)
    if size(QR_input, 2) ~= 75
        error('Input dimension mismatch: expected N×75 features (Q+R).');
    end
    y_pred = predict(rfModel, QR_input);
    fprintf('Predicted %d mask values.\n', numel(y_pred));
end
