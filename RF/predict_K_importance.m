function K_importance = predict_K_importance(A, B, adjMtx, modelPath)
% predict_K_importance - Predicts K element-wise importance using trained RF
%
% Inputs:
%   A        - System state matrix (Nx × Nx)
%   B        - Input matrix (Nx × Nu)
%   adjMtx   - Adjacency matrix of the graph (Nnode × Nnode)
%   modelPath - Full path to the trained TreeBagger .mat file
%
% Output:
%   K_importance - Predicted importance matrix (Nu × Nx)

    % --- Load model ---
    data = load(modelPath);
    rfModel = data.rfModel;

    % --- Dimensions ---
    [Nx, Nu] = size(B);

    % --- Feature Construction ---
    feat_d      = reshape(adjMtx, 1, []);  % vectorized distances
    feat_colA   = sum(abs(A), 1);          % column norm of A
    feat_diagA  = diag(A)';                % A diagonals
    feat_bPow   = sum(B.^2, 1);            % actuator strength
    feat_align  = dot(B, A, 1);            % alignment feature

    % --- Repeat and align features for all (i,j) pairs ---
    feat_mat = [];
    for i = 1:Nu
        for j = 1:Nx
            dij = adjMtx(i,j);                  % distance between actuator i and state j
            feat_entry = [dij, feat_colA(j), feat_diagA(j), feat_bPow(i), feat_align(i)];
            feat_mat = [feat_mat; feat_entry];  % one row per K(i,j)
        end
    end

    % --- Predict using RF model ---
    pred_vec = predict(rfModel, feat_mat);  % (Nu*Nx × 1)
    K_importance = reshape(pred_vec, Nu, Nx);

end
