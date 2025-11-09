function print_top_list(currGen, bestGens, bestCosts, bestChromsCts)

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
