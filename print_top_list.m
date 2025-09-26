function print_top_list(currGen, bestGens, bestCosts, bestChromsCts)
% PRINT_TOP_LIST 打印当前种群的最优个体列表信息
%
% 输入参数：
%   currGen        — 当前代数（整数）
%   bestGens       — 长度为 nTop 的向量，记录每个最优个体所在的代数
%   bestCosts      — 长度为 nTop 的向量，记录对应个体的成本
%   bestChromsCts  — 长度为 nTop 的 cell 数组，每个元素是 [# sensors, # actuators]

    fprintf('\n==== Generation %d ====\n', currGen);
    fprintf('%4s | %10s | %10s | %10s\n', 'Rank', 'Gen', 'Cost', '[Sensors, Actuators]');
    fprintf('----------------------------------------------\n');
    for i = 1:length(bestCosts)
        cts = bestChromsCts{i};
        fprintf('%4d | %10d | %10.4f | [%2d, %2d]\n', ...
                i, bestGens(i), bestCosts(i), cts(1), cts(2));
    end
    fprintf('==============================================\n');
end
