
function mx = mutate_cts(x, pMutate, Nx, Nu)
    mx = x;
    p = rand(1,2);
    if p(1) < pMutate
        mx(1) = randi(Nx);
    end
    if p(2) < pMutate
        mx(2) = randi(Nu);
    end
end