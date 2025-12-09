function Q = critic_forward(critic, S, Z)
    % S: sDim × batch
    % Z: dAct × batch
    X = [S; Z];

    h1 = critic.W1 * X + critic.b1;
    a1 = max(h1,0);

    h2 = critic.W2 * a1 + critic.b2;
    a2 = max(h2,0);

    Q  = critic.W3 * a2 + critic.b3;   % 1 × batch
end
