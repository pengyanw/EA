function [child1, child2] = reproduce_discrete(parent1, parent2)
% point intersection
    n = length(parent1);
    crossoverPoint = randi(n);
    child1 = [parent1(1:crossoverPoint), parent2(crossoverPoint+1:end)];
    child2 = [parent2(1:crossoverPoint), parent1(crossoverPoint+1:end)];
    %printf('child1:',child1)
end