function actor = init_actor(sDim, dAct)

    actor.W1 = randn(128, sDim) * 0.1;
    actor.b1 = zeros(128, 1);

    actor.W2 = randn(128, 128) * 0.1;
    actor.b2 = zeros(128, 1);

    actor.W3 = randn(dAct, 128) * 0.1;
    actor.b3 = zeros(dAct, 1);

    % Adam 状态
    actor.mW1 = zeros(size(actor.W1));
    actor.mb1 = zeros(size(actor.b1));
    actor.mW2 = zeros(size(actor.W2));
    actor.mb2 = zeros(size(actor.b2));
    actor.mW3 = zeros(size(actor.W3));
    actor.mb3 = zeros(size(actor.b3));

    actor.vW1 = zeros(size(actor.W1));
    actor.vb1 = zeros(size(actor.b1));
    actor.vW2 = zeros(size(actor.W2));
    actor.vb2 = zeros(size(actor.b2));
    actor.vW3 = zeros(size(actor.W3));
    actor.vb3 = zeros(size(actor.b3));
end

