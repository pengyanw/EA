function critic = init_critic(sDim, dAct)

    critic.W1 = randn(128, sDim + dAct) * 0.1;
    critic.b1 = zeros(128,1);

    critic.W2 = randn(128,128) * 0.1;
    critic.b2 = zeros(128,1);

    critic.W3 = randn(1,128) * 0.1;
    critic.b3 = 0;

    % Adam 状态
    critic.mW1 = zeros(size(critic.W1));
    critic.mb1 = zeros(size(critic.b1));
    critic.mW2 = zeros(size(critic.W2));
    critic.mb2 = zeros(size(critic.b2));
    critic.mW3 = zeros(size(critic.W3));
    critic.mb3 = zeros(size(critic.b3));

    critic.vW1 = zeros(size(critic.W1));
    critic.vb1 = zeros(size(critic.b1));
    critic.vW2 = zeros(size(critic.W2));
    critic.vb2 = zeros(size(critic.b2));
    critic.vW3 = zeros(size(critic.W3));
    critic.vb3 = zeros(size(critic.b3));
end
