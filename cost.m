function [costs,currK] = cost(currPopBin, currPopCts, A, B_, C_, Q, R_, K_, baseCost, sensPen, actPen, commPen)
    [YIdx, UIdx] = postprocess(currPopBin, currPopCts, K_);

    popSize = length(currPopBin);
    costs   = nan(popSize, 1);
    currK = [];
    for i = 1:popSize
        KSupp     = reshape(currPopBin{i}, size(K_));
        K         = K_;
        K(~KSupp) = 0;
        K(:, ~YIdx{i}) = []; % Eliminate unsensed columns
        K(~UIdx{i}, :) = []; % Eliminate unactuated rows

        K1         = K_;
        K1(~KSupp) = 0;
        K1(:, ~YIdx{i}) = 0; % Eliminate unsensed columns
        K1(~UIdx{i}, :) = 0; % Eliminate unactuated rows

        currK(i,:,:)  = K1;
        B = B_(:, UIdx{i});
        C = C_(YIdx{i}, :);
        R = R_(UIdx{i}, UIdx{i});
        
        lqrCost  = get_lqr_cost(A, B, Q, R, K*C) ./ baseCost; % normalized

        costs(i) = lqrCost + sensPen*currPopCts{i}(1) + actPen*currPopCts{i}(2) + commPen*nnz(K);        
    end
end