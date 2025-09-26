function mutant = mutate_discrete(x, pMutate)
% 以 pMutate 概率对二进制数组每一位进行翻转
    mutant = x;
    for i = 1:length(x)
        if rand < pMutate
            mutant(i) = ~x(i); % 0 <-> 1
        end
    end
end