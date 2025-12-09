function grad = critic_backward(critic, sBatch, zBatch, diffQ)
    % diffQ: (qPred - target) of size [1, batch]
    batch = size(sBatch,2);

    grad.W1 = zeros(size(critic.W1));
    grad.b1 = zeros(size(critic.b1));
    grad.W2 = zeros(size(critic.W2));
    grad.b2 = zeros(size(critic.b2));
    grad.W3 = zeros(size(critic.W3));
    grad.b3 = 0;

    for t = 1:batch
        s = sBatch(:,t);
        z = zBatch(:,t);
        x0 = [s; z];

        % forward (store intermediates)
        h1 = critic.W1*x0 + critic.b1; a1 = max(h1,0);
        h2 = critic.W2*a1 + critic.b2; a2 = max(h2,0);
        q  = critic.W3*a2 + critic.b3;

        dq = diffQ(t);

        % layer3
        grad.W3 = grad.W3 + dq * a2';
        grad.b3 = grad.b3 + dq;

        da2 = critic.W3' * dq;
        dh2 = da2 .* (h2>0);

        % layer2
        grad.W2 = grad.W2 + dh2 * a1';
        grad.b2 = grad.b2 + dh2;

        da1 = critic.W2' * dh2;
        dh1 = da1 .* (h1>0);

        % layer1
        grad.W1 = grad.W1 + dh1 * x0';
        grad.b1 = grad.b1 + dh1;
    end

    % average
    f = 1/batch;
    grad.W1 = grad.W1 * f;
    grad.b1 = grad.b1 * f;
    grad.W2 = grad.W2 * f;
    grad.b2 = grad.b2 * f;
    grad.W3 = grad.W3 * f;
    grad.b3 = grad.b3 * f;
end
