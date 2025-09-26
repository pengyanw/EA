function repProb = generate_reproduction_probs(costs, lqrcost)
    % 输入:
    %   costs : 个体的总 cost (Inf 或 大于阈值表示不稳定)
    % 输出:
    %   repProb : reproduction probability 向量

    % 超参数
    penalty = 1e-6;  % 不稳定个体的最小权重
    n = length(costs);
    repProb = zeros(n,1);

    % 判断稳定性 (这里假定 cost >= 1e6 或 Inf 代表不稳定)
    stableMask   = isfinite(costs) & (log(costs/lqrcost) < 5);
    unstableMask = ~stableMask;

    % ===== 稳定个体 =====
    if any(stableMask)
        % 转换为越小越好的“适应度” (cost 越小适应度越大)
        fit = max(costs(stableMask)) - costs(stableMask) + 1e-6;
        repProb(stableMask) = fit;
    end

    % ===== 不稳定个体 =====
    repProb(unstableMask) = penalty;

    % ===== 归一化 =====
    if sum(repProb) == 0
        repProb(:) = 1/n;  % fallback: 平均分配
    else
        repProb = repProb / sum(repProb);
    end
end
