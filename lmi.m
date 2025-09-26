% 示例系统（你可替换为真实系统）
A = [1.0 0; 0 1.1];   % 不稳定
B = [1.0; 0.5];       % 单输入
Q = eye(2);
R = 1;

% 基础 LQR 解作为结构参考
K_lqr = -dlqr(A, B, Q, R);

% bin 结构：根据 LQR 的非零结构初始化
[n, m] = size(B);
bin = reshape(abs(K_lqr) > 1e-5, 1, []);

fprintf('=== Testing feasible cts for given bin ===\n');
max_r = n;  % 最大行
max_c = m;  % 最大列

success_cnt = 0;

for r = 0:max_r
    for c = 0:max_c
        cts = [r, c];
        K_try = lmi_repair(A, B, Q, R, bin, cts);
        if ~isempty(K_try)
            fprintf('✅ Feasible for cts = [%d, %d]\n', r, c);
            success_cnt = success_cnt + 1;
        else
            fprintf('❌ Infeasible for cts = [%d, %d]\n', r, c);
        end
    end
end

fprintf('Total feasible (bin fixed): %d combinations\n', success_cnt);

