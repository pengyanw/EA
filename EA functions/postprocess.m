% Updated Postprocess using L1-norm strength (not just count)
function [YIdx, UIdx] = postprocess(currPopBin, currPopCts, K_)
    popSize  = length(currPopBin);
    YIdx     = cell(popSize, 1);
    UIdx     = cell(popSize, 1);
    [Nu, Nx] = size(K_);

    for i = 1:popSize
        KSupp = reshape(currPopBin{i}, size(K_));
        K_masked = abs(K_ .* KSupp);  % masked K, take absolute for scoring

        % ==== Sensor (Column) Importance ====
        col_strength = sum(K_masked.^2, 1);  % 1×Nx
        YIdx{i} = false(1, Nx);
        [~, sortedCols] = sort(col_strength, 'descend');
        topCols = sortedCols(1:min(currPopCts{i}(1), Nx));
        YIdx{i}(topCols) = true;

        % ==== Actuator (Row) Importance ====
        row_strength = sum(K_masked.^2, 2);  % Nu×1
        UIdx{i} = false(1, Nu);
        [~, sortedRows] = sort(row_strength, 'descend');
        topRows = sortedRows(1:min(currPopCts{i}(2), Nu));
        UIdx{i}(topRows) = true;
    end
end
