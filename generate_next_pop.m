function [nextPopBin, nextPopCts] = generate_next_pop(currPopBin, currPopCts, repProb, pMutate, Nx, Nu)
    popSize    = length(currPopBin);
    nextPopBin = cell(popSize, 1);
    nextPopCts = cell(popSize, 1);

    for j=1:popSize/2
        is = randsample(1:popSize, 2, true, repProb);
        while is(1) == is(2) % Do not allow self-reproduction
            is = randsample(1:popSize, 2, true, repProb);
        end
        % Generate the discrete-valued parts of next generation
        [c1, c2] = reproduce_discrete(currPopBin{is(1)}, currPopBin{is(2)});
        nextPopBin{2*j-1} = mutate_discrete(c1, pMutate);
        nextPopBin{2*j}   = mutate_discrete(c2, pMutate);

        % Generate the continuous-valued parts of next generation
        [c1, c2] = reproduce_continuous(currPopCts{is(1)}, currPopCts{is(2)});
        nextPopCts{2*j-1} = round(mutate_cts(c1, pMutate, Nx, Nu));
        nextPopCts{2*j}   = round(mutate_cts(c2, pMutate, Nx, Nu));
    end
end
% default function
function [child1, child2] = reproduce_discrete(parent1, parent2)
% point intersection
    n = length(parent1);
    crossoverPoint = randi(n);
    child1 = [parent1(1:crossoverPoint), parent2(crossoverPoint+1:end)];
    child2 = [parent2(1:crossoverPoint), parent1(crossoverPoint+1:end)];
end