function [child1, child2] = reproduce_continuous(parent1, parent2)

% parent2 = max(parent1, parent2);
% parent1 = min(parent1, parent2);
    alpha = rand;
    child1 = round(alpha * parent1 + (1 - alpha) * parent2);
    child2 = round((1 - alpha) * parent1 + alpha * parent2);
end
