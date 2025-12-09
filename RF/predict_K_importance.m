function K_map = predict_K_importance(A, B, adjMtx, modelFile)
% Predict |K| importance map for a new system using trained RF model

% Load RF model
loaded = load(modelFile, 'rfModel');
rfModel = loaded.rfModel;

% Dimensions
[Nx, Nu] = size(B);
numNodes = size(adjMtx, 1);
assert(mod(Nx,numNodes)==0 && mod(Nu,numNodes)==0, 'Nx/Nu not divisible by #nodes');

nodePerAct  = Nu / numNodes;
nodePerState = Nx / numNodes;

% --- Distance Matrix Expansion ---
G = graph(adjMtx);
D_node = distances(G);
D_node(isinf(D_node)) = max(D_node(~isinf(D_node))) + 1;
D_exp = kron(D_node, ones(nodePerAct, nodePerState));  % Nu × Nx

% --- Global norms ---
normA = norm(A, 'fro');
normB = norm(B, 'fro');

% --- Local AB features ---
colA  = sum(abs(A), 1)';      % Nx × 1
diagA = abs(diag(A));         % Nx × 1
bPow  = vecnorm(B, 2, 1)';    % Nu × 1
Galign = abs(B' * A);         % Nu × Nx

% --- Assemble features ---
[II, JJ] = ndgrid(1:Nu, 1:Nx);  % actuator-state pairs
feat_d      = D_exp(:);
feat_colA   = colA(JJ(:));
feat_diagA  = diagA(JJ(:));
feat_bPow   = bPow(II(:));
feat_align  = Galign(:);

X_pred = [feat_d, feat_colA, feat_diagA, feat_bPow, feat_align];

% --- Predict
Y_pred = predict(rfModel, X_pred);    % length = Nu*Nx
K_map = reshape(Y_pred, Nu, Nx);      % Nu × Nx

end
