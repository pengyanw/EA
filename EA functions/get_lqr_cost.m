function val = get_lqr_cost(A, B, Q, R, K)
    A_ = (A+B*K)';
    Q_ = Q + K'*R*K;
    
    if ~isempty(find(abs(eig(A_)) >= 1, 1))
        val = inf; 
        return;
    end

    try
        P = dlyap(A_, Q_);
    catch
        val = inf;
        return;
    end

    val = norm(P, 'fro').^2;
end