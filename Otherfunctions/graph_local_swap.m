function K_new = graph_local_swap(A, B, Q, R, Adj, K, num_links, D)
%--------------------------------------------------------------------------
% Graph-Guided Local Structural Surgery (GG-LSS)
%
%  输入:
%     A,B      : 系统矩阵
%     Q,R      : LQR 权重
%     Adj      : adjacency matrix (Nx x Nx)
%     K        : 当前 sparse K
%     num_links: 当前 sparsity 预算
%     D        : (可选) shortest-path matrix
%
%  输出:
%     K_new    : 修改后 K（若不稳定自动 fallback）
%
%  说明:
%   - 本函数仅修改 K 的一行，实现 “删远加近” 的拓扑结构跳变
%   - 核心操作: 删除 row i 上最远的非零项，添加邻居中最近的零项
%   - 增强稳定性: 用 row-wise top-k + kernel gating 重建本行
%--------------------------------------------------------------------------

Nu = size(B,2);
Nx = size(A,1);

% 如果没有 shortest-path matrix，则计算
if nargin < 8 || isempty(D)
    G = graph(Adj);
    D = distances(G);
end

%==== Step 0: 计算 dense LQR 作为基准方向 ================================
try
    K_dense = -dlqr(A, B, Q, R);
catch
    % 如果 LQR 不可解
    K_new = K;
    return;
end

%==== Step 1: 选择一行 i（重点：偏向 diagonal 强的行） ====================
diag_strength = abs(diag(K_dense(1:min(Nu,Nx), 1:min(Nu,Nx))));
prob = diag_strength / (sum(diag_strength) + eps);
i = randsample(1:Nu,1,true,prob);

row = K(i,:);
nz_idx = find(row ~= 0);

if isempty(nz_idx)
    K_new = K; 
    return;
end

%==== Step 2: 删除 row i 中 “最远的连接” ================================
[~, id_far] = max(D(i, nz_idx));
j_far = nz_idx(id_far);

K2 = K;
K2(i, j_far) = 0;

%==== Step 3: 添加 row i 中最邻近且当前为 0 的连接 =======================
zero_idx = find(row == 0);

if isempty(zero_idx)
    K_new = K;
    return;
end

[~, id_near] = min(D(i, zero_idx));
j_near = zero_idx(id_near);

% 添加基于 dense K 的方向
K2(i, j_near) = K_dense(i, j_near);

%==== Step 4: 用 kernel gating + row-wise top-k 修理本行 ===================
% kernel: exponential decay by graph distance
beta = 0.7;
W = exp(-beta * D(i,:));

row_weighted = K_dense(i,:) .* W;

% 只修 row i
[~, idx_sort] = sort(abs(row_weighted), 'descend');
k_row = min(num_links, Nx);
keep = idx_sort(1:k_row);

new_row = zeros(1,Nx);
new_row(keep) = row_weighted(keep);

K2(i,:) = new_row;

% 归一化 row（保持能量一致）
nr = norm(K2(i,:),2) + 1e-12;
nd = norm(K_dense(i,:),2) + 1e-12;
K2(i,:) = K2(i,:) * (nd / nr);

%==== Step 5: 稳定性检查 + fallback ======================================
Acl = A + B*K2;
if max(abs(eig(Acl))) < 1
    K_new = K2;
    return;
end

% fallback：线搜索 η
eta = 1.0;
K_new = K;
while eta > 1e-3
    K_mix = (1-eta)*K + eta*K2;
    if max(abs(eig(A + B*K_mix))) < 1
        K_new = K_mix;
        return;
    end
    eta = eta / 2;
end

% fallback failed → 返回原 K
K_new = K;

end
