function repProb = generate_reproduction_probs_origin(costs, lqrcosts)
    % Weight, adjust sensitivity to LQR
    alpha = 1.0;   % weight for total cost
    beta  = 10.0;   % weight for LQR cost (higher = more sensitive)

    
    cost_fit = max(costs) - costs + 1e-6;
    lqr_fit  = max(lqrcosts) - lqrcosts + 1e-6;

    % integrate fitness
    fitness = alpha * cost_fit + beta * lqr_fit;
    fitness = max(fitness, 0);  % nonnegative

    % transform to reproduction probabilities
    if all(fitness == 0)
        repProb = ones(size(costs)) / length(costs);  % fallback
    else
        repProb = fitness / sum(fitness);
    end
end
