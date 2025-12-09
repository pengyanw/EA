function loss = lyapunov_loss(dlnetV, dlnetU, Xbatch)
    alpha = 0.1;

    V = forward(dlnetV, Xbatch);
    U = forward(dlnetU, Xbatch);

    X = extractdata(Xbatch);
    Uu = extractdata(U);

    n = size(X,2);
    f = zeros(2,n);
    for k = 1:n
        f(:,k) = pendulum_dyn(X(:,k), Uu(k));
    end
    f = dlarray(f,'CB');

    % gradient wrt Xbatch
    dVdx = dlgradient(sum(V), Xbatch);

    dVdt = sum(dVdx .* f, 1);

    penalty = relu(dVdt + alpha * sum(Xbatch.^2,1));
    loss = mean(penalty);
end
